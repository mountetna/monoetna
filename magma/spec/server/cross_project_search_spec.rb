describe 'cross project search routes' do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  before do
    %w[subject biospecimen sc_seq sc_seq_pool project assay other_assay].each { |name| add_model('coprojects_template', name) }
  end

  after do
    (@loaded_projects || []).each { |project_name| unload_project(project_name) }
  end

  # Inserts the model row, or updates its template mapping when it already exists.
  def add_model(project_name, model_name, template_model_name: nil)
    values = {
      template_project_name: template_model_name && 'coprojects_template',
      template_model_name: template_model_name
    }
    dataset = Magma.instance.db[:models].where(project_name: project_name, model_name: model_name)
    return dataset.update(values) if dataset.first

    Magma.instance.db[:models].insert(values.merge(project_name: project_name, model_name: model_name))
  end

  # Inserts the attribute row, or updates it when it already exists.
  def create_attribute(project_name, model_name, attribute_name, type: 'string', link_model_name: nil, column_name: attribute_name, link_attribute_name: nil)
    values = {
      type: type,
      column_name: column_name,
      link_model_name: link_model_name,
      link_attribute_name: link_attribute_name
    }
    dataset = Magma.instance.db[:attributes].where(
      project_name: project_name,
      model_name: model_name,
      attribute_name: attribute_name
    )
    return dataset.update(values) if dataset.first

    Magma.instance.db[:attributes].insert(values.merge(
      project_name: project_name,
      model_name: model_name,
      attribute_name: attribute_name
    ))
  end

  # Replaces the model's table with one holding the given columns and rows.
  def add_table(project_name, model_name, columns, rows = [])
    table = Sequel[project_name.to_sym][model_name.pluralize.to_sym]
    db = Magma.instance.db
    db.create_schema(project_name.to_sym, if_not_exists: true)
    db.drop_table?(table)
    db.create_table(table) do
      primary_key :id
      columns.each do |name, type|
        column name, type == :string ? String : type
      end
    end
    rows.each { |row| db[table].insert(row) }
  end

  def add_parent(project_name, model_name, parent_model)
    create_attribute(
      project_name,
      model_name,
      parent_model,
      type: 'parent',
      link_model_name: parent_model,
      column_name: "#{parent_model}_id",
      link_attribute_name: model_name
    )
  end

  def add_collection(project_name, model_name, child_model, parent_attribute)
    create_attribute(
      project_name,
      model_name,
      child_model,
      type: 'collection',
      link_model_name: child_model,
      column_name: child_model,
      link_attribute_name: parent_attribute
    )
  end

  def add_search_lineage(project_name, intermediate: nil)
    add_model(project_name, 'patient', template_model_name: 'subject')
    add_model(project_name, 'specimen', template_model_name: 'biospecimen')
    add_model(project_name, 'tube', template_model_name: 'sc_seq')
    add_model(project_name, 'pool', template_model_name: 'sc_seq_pool')
    add_parent(project_name, 'tube', 'specimen')
    create_attribute(project_name, 'tube', 'pool', type: 'link', link_model_name: 'pool', column_name: 'pool_id')

    if intermediate
      add_model(project_name, intermediate)
      add_parent(project_name, 'specimen', intermediate)
      add_collection(project_name, 'specimen', 'tube', 'specimen')
      add_parent(project_name, intermediate, 'patient')
      add_collection(project_name, intermediate, 'specimen', intermediate)
      add_collection(project_name, 'patient', intermediate, 'patient')
      parent_column = "#{intermediate}_id".to_sym
      add_table(project_name, 'patient', {}, [{}])
      add_table(project_name, intermediate, { patient_id: :integer }, [{ patient_id: 1 }])
      add_table(project_name, 'specimen', { parent_column => :integer }, [{ parent_column => 1 }])
    else
      add_parent(project_name, 'specimen', 'patient')
      add_collection(project_name, 'specimen', 'tube', 'specimen')
      add_collection(project_name, 'patient', 'specimen', 'patient')
      add_table(project_name, 'patient', {}, [{}])
      add_table(project_name, 'specimen', { patient_id: :integer }, [{ patient_id: 1 }])
    end

    add_table(project_name, 'pool', {}, [{}])
    add_table(project_name, 'tube', { specimen_id: :integer, pool_id: :integer }, [{ specimen_id: 1, pool_id: 1 }])
  end

  # Projects created by raw inserts are not in Magma's project cache; load
  # them fresh and drop them again so no stale model classes leak between examples.
  def load_projects(*project_names)
    @loaded_projects = project_names
    project_names.each do |project_name|
      unload_project(project_name)
      Magma.instance.get_or_load_project(project_name)
    end
  end

  def unload_project(project_name)
    Magma.instance.magma_projects.delete(project_name.to_sym)
    constant = project_name.to_s.camel_case.to_sym
    Object.send(:remove_const, constant) if Object.const_defined?(constant, false)
  end

  def superviewer_header
    header(*Etna::TestAuth.token_header(
      email: 'argus@twelve-labors.org', name: 'Argus', perm: 'v:administration', exp: Time.now.to_i + 6000
    ))
  end

  def post_as_superuser(endpoint, body = {})
    auth_header(:superuser)
    json_post(endpoint, body)
    last_response
  end

  def add_field(label, query, position: nil, kind: nil, modality_id: :none)
    json = { label: label, query: query }
    json[:modality_id] = modality_id unless modality_id == :none
    json[:kind] = kind || ((Array(query).flatten & %w[::distinct ::first]).any? ? 'filter' : 'count')
    json[:position] = position if position
    post_as_superuser('cross_project_search/fields', json)
  end

  def update_field(id, body)
    post_as_superuser("cross_project_search/fields/#{id}", body)
  end

  def delete_field(id)
    post_as_superuser("cross_project_search/fields/#{id}/delete")
  end

  def add_modality(name, model = 'assay', **extra)
    post_as_superuser('cross_project_search/modalities', { name: name, model: model }.merge(extra))
  end

  def update_modality(id, body)
    post_as_superuser("cross_project_search/modalities/#{id}", body)
  end

  def delete_modality(id)
    post_as_superuser("cross_project_search/modalities/#{id}/delete")
  end

  def fields
    auth_header(:superuser)
    get('/cross_project_search/fields')
    expect(last_response.status).to eq(200), last_response.body
    json_body[:fields]
  end

  def modalities
    auth_header(:superuser)
    get('/cross_project_search/modalities')
    expect(last_response.status).to eq(200), last_response.body
    json_body[:modalities]
  end

  def run_queries(body)
    post_as_superuser('cross_project_search/query', body)
  end

  def results(body)
    response = run_queries(body)
    expect(response.status).to eq(200), response.body
    json_body[:projects]
  end

  def error_message
    json_body[:error]
  end

  context 'access' do
    it 'rejects unauthenticated requests' do
      get('/cross_project_search/fields')
      expect(last_response.status).to eq(401)
      get('/cross_project_search/modalities')
      expect(last_response.status).to eq(401)
      json_post('cross_project_search/fields', label: 'x', query: ['subject', '::count'])
      expect(last_response.status).to eq(401)
      json_post('cross_project_search/modalities', name: 'x', model: 'assay')
      expect(last_response.status).to eq(401)
      json_post('cross_project_search/query', queries: [['subject', '::count']])
      expect(last_response.status).to eq(401)
    end

    it 'rejects a project viewer on every route' do
      auth_header(:viewer)
      get('/cross_project_search/fields')
      expect(last_response.status).to eq(403)
      get('/cross_project_search/modalities')
      expect(last_response.status).to eq(403)
      json_post('cross_project_search/fields', label: 'x', query: ['subject', '::count'])
      expect(last_response.status).to eq(403)
      json_post('cross_project_search/modalities', name: 'x', model: 'assay')
      expect(last_response.status).to eq(403)
      json_post('cross_project_search/query', queries: [['subject', '::count']])
      expect(last_response.status).to eq(403)
      json_post('cross_project_search/fields/1/delete', {})
      expect(last_response.status).to eq(403)
      json_post('cross_project_search/modalities/1/delete', {})
      expect(last_response.status).to eq(403)
    end

    it 'lets a superviewer read and search but not write' do
      superviewer_header
      get('/cross_project_search/fields')
      expect(last_response.status).to eq(200)
      get('/cross_project_search/modalities')
      expect(last_response.status).to eq(200)
      json_post('cross_project_search/query', queries: [['subject', '::count']])
      expect(last_response.status).to eq(200), last_response.body

      json_post('cross_project_search/fields', label: 'x', query: ['subject', '::count'])
      expect(last_response.status).to eq(403)
      json_post('cross_project_search/fields/1', label: 'x')
      expect(last_response.status).to eq(403)
      json_post('cross_project_search/fields/1/delete', {})
      expect(last_response.status).to eq(403)
      json_post('cross_project_search/modalities', name: 'x', model: 'assay')
      expect(last_response.status).to eq(403)
      json_post('cross_project_search/modalities/1', name: 'x')
      expect(last_response.status).to eq(403)
      json_post('cross_project_search/modalities/1/delete', {})
      expect(last_response.status).to eq(403)
    end
  end

  describe SearchFieldsController do
    it 'starts empty and creates, lists, updates and deletes rows as a superuser' do
      expect(fields).to eq([])

      response = add_field('Subjects', ['subject', '::count'])
      expect(response.status).to eq(200), response.body
      created = json_body[:field]
      expect(created).to include(label: 'Subjects', query: ['subject', '::count'], kind: 'count', modality_id: nil, modality: nil, position: 1)
      expect(created[:updated_at]).not_to be_nil
      expect(fields.map { |field| field[:id] }).to eq([created[:id]])

      update_field(created[:id], label: 'Subject count')
      expect(last_response.status).to eq(200), last_response.body
      expect(json_body[:field]).to include(label: 'Subject count', query: ['subject', '::count'])

      update_field(created[:id], kind: 'filter', query: ['subject', '::distinct', 'species'])
      expect(json_body[:field]).to include(kind: 'filter', query: ['subject', '::distinct', 'species'])

      expect(delete_field(created[:id]).status).to eq(200)
      expect(fields).to eq([])
    end

    it 'stores the query as json' do
      add_field('Samples', ['biospecimen', ['::has', 'x'], '::count'])
      expect(Magma.instance.db.schema(:search_fields).to_h[:query][:db_type]).to eq('json')
      expect(fields.first[:query]).to eq(['biospecimen', ['::has', 'x'], '::count'])
    end

    it 'changes updated_at on edit and leaves created_at' do
      add_field('Subjects', ['subject', '::count'])
      id = json_body[:field][:id]
      Magma.instance.db[:search_fields].where(id: id).update(updated_at: Time.now - 3600)
      before_edit = fields.first[:updated_at]
      created_at = Magma.instance.db[:search_fields].where(id: id).get(:created_at)

      update_field(id, label: 'Renamed')
      expect(Time.parse(fields.first[:updated_at])).to be > Time.parse(before_edit)
      expect(Magma.instance.db[:search_fields].where(id: id).get(:created_at)).to eq(created_at)
    end

    it 'puts new rows at the end and orders by position then id' do
      add_field('A', ['subject', '::count'])
      add_field('B', ['subject', '::count'])
      add_field('C', ['subject', '::count'], position: 0)
      expect(fields.map { |field| field[:label] }).to eq(%w[C A B])

      update_field(fields.first[:id], position: 99)
      expect(fields.map { |field| field[:label] }).to eq(%w[A B C])
    end

    context 'modality_id' do
      it 'is null or an existing modality, and the list names it' do
        add_field('Subjects', ['subject', '::count'])
        expect(json_body[:field]).to include(modality_id: nil, modality: nil)

        expect(add_field('A', ['subject', '::count'], modality_id: 9999).status).to eq(422)
        expect(error_message).to match(/does not exist/)
        expect(add_field('A', ['subject', '::count'], modality_id: 'single cell').status).to eq(422)
        expect(error_message).to match(/integer/)

        add_modality('single cell', 'assay')
        modality_id = json_body[:modality][:id]
        add_field('Kind', ['subject', '::distinct', 'species'], modality_id: modality_id)
        expect(json_body[:field]).to include(kind: 'filter', modality_id: modality_id, modality: 'single cell')

        id = json_body[:field][:id]
        update_field(id, modality_id: '')
        expect(json_body[:field]).to include(modality_id: nil, modality: nil)
        update_field(id, modality_id: modality_id)
        expect(json_body[:field]).to include(modality_id: modality_id, modality: 'single cell')
        expect(fields.map { |field| field[:modality] }).to eq([nil, 'single cell'])
      end
    end

    context 'kind' do
      it 'round-trips an explicit kind and lists it' do
        add_field('Subjects', ['subject', '::count'], kind: 'count')
        expect(json_body[:field]).to include(kind: 'count')
        add_field('Species', ['subject', '::distinct', 'species'], kind: 'filter')
        expect(json_body[:field]).to include(kind: 'filter')
        expect(fields.map { |field| field[:kind] }).to eq(%w[count filter])
      end

      it 'returns 422 when the kind and the verb disagree' do
        shared = /A count must end in ::count\. A filter must end in ::distinct or ::first, followed by an attribute\./
        {
          'count without ::count' => ['subject', '::distinct', 'species'],
          'filter with the wrong verb' => ['subject', '::count'],
          'filter without an attribute' => ['subject', '::first']
        }.each do |name, query|
          kind = name.start_with?('count') ? 'count' : 'filter'
          expect(add_field(name, query, kind: kind).status).to eq(422)
          expect(error_message).to match(shared)
        end
        expect(fields).to eq([])
      end

      it 'keeps the kind in the table to count or filter' do
        add_field('A', ['subject', '::count'])
        expect do
          Magma.instance.db.transaction(savepoint: true) do
            Magma.instance.db[:search_fields].where(id: json_body[:field][:id]).update(kind: 'column')
          end
        end.to raise_error(Sequel::CheckConstraintViolation)
      end

      it 'checks kind against the query on update' do
        add_field('A', ['subject', '::count'])
        id = json_body[:field][:id]
        expect(update_field(id, query: ['subject', '::distinct', 'species']).status).to eq(422)
        expect(update_field(id, kind: 'filter').status).to eq(422)
        expect(update_field(id, label: 'Still fine').status).to eq(200)
        expect(json_body[:field][:kind]).to eq('count')
      end
    end

    it 'returns 404 when updating or deleting a missing field' do
      expect(update_field(9999, label: 'x').status).to eq(404)
      expect(error_message).to eq('No such field')
      expect(delete_field(9999).status).to eq(404)
    end

    it 'stores a query whose grammar Magma judges when it runs' do
      expect(add_field('Unmapped', ['biospecimen', ['::has', 'x'], '::count']).status).to eq(200)
      expect(add_field('Nested', ['subject', ['a', '::bogus'], '::count']).status).to eq(200)
      expect(add_field('Missing', ['nope', '::count']).status).to eq(200)
      expect(add_field('Extra', ['subject', '::count', 'species']).status).to eq(200)
    end

    {
      'an empty label' => [{ label: ' ', query: ['subject', '::count'] }, /label/i],
      'a missing label' => [{ query: ['subject', '::count'] }, /label/i],
      'a query that is not an array' => [{ label: 'x', query: 'subject' }, /non-empty JSON array/],
      'a query that is an object' => [{ label: 'x', query: { a: 1 } }, /non-empty JSON array/],
      'an empty array' => [{ label: 'x', query: [] }, /non-empty JSON array/],
      'a query that does not start with a name' => [{ label: 'x', query: [1, '::count'] }, /non-empty JSON array/],
      'no result verb' => [{ label: 'x', kind: 'count', query: ['subject', 'species'] }, /A count must end in ::count/],
      'a bad position' => [{ label: 'x', kind: 'count', query: ['subject', '::count'], position: 'two' }, /integer/],
      'a missing kind' => [{ label: 'x', query: ['subject', '::count'] }, /kind/i],
      'an unknown kind' => [{ label: 'x', kind: 'column', query: ['subject', '::count'] }, /'count' or 'filter'/]
    }.each do |name, (body, message)|
      it "rejects #{name} with a 422 message on create" do
        expect(post_as_superuser('cross_project_search/fields', body).status).to eq(422), last_response.body
        expect(error_message).to match(message)
        expect(fields).to eq([])
      end
    end

    it 'rejects bad edits with a 422 and leaves the row unchanged' do
      add_field('Subjects', ['subject', '::count'])
      id = json_body[:field][:id]

      expect(update_field(id, query: 'subject').status).to eq(422)
      expect(error_message).to match(/non-empty JSON array/)
      expect(update_field(id, label: '').status).to eq(422)
      expect(fields.first).to include(label: 'Subjects', query: ['subject', '::count'])
    end
  end

  describe SearchModalitiesController do
    it 'creates, lists, updates and deletes modalities' do
      expect(modalities).to eq([])

      response = add_modality('single cell', 'assay')
      expect(response.status).to eq(200), response.body
      created = json_body[:modality]
      expect(created).to include(name: 'single cell', model: 'assay', position: 1)
      expect(created).not_to have_key(:query)
      expect(created[:updated_at]).not_to be_nil

      update_modality(created[:id], name: 'Single cell', model: 'other_assay')
      expect(last_response.status).to eq(200), last_response.body
      expect(json_body[:modality]).to include(name: 'Single cell', model: 'other_assay')
      expect(modalities.map { |modality| modality[:id] }).to eq([created[:id]])

      expect(delete_modality(created[:id]).status).to eq(200)
      expect(modalities).to eq([])
    end

    it 'requires the model to be a template model' do
      expect(add_modality('flow', 'not_a_real_model').status).to eq(422)
      expect(error_message).to match(/Unknown model 'not_a_real_model'/)
      add_modality('flow', 'assay')
      update_modality(json_body[:modality][:id], model: 'nope')
      expect(last_response.status).to eq(422)
      expect(modalities.map { |modality| modality[:model] }).to eq(['assay'])
    end

    it 'requires a name and a non-blank model' do
      expect(add_modality('  ', 'assay').status).to eq(422)
      expect(error_message).to match(/name/i)
      expect(add_modality('flow', '   ').status).to eq(422)
      expect(error_message).to match(/model/i)
      expect(add_modality('flow', nil).status).to eq(422)
      expect(error_message).to match(/model/i)
      expect(modalities).to eq([])
    end

    it 'keeps names unique, ignoring case' do
      add_modality('single cell', 'assay')
      expect(add_modality(' Single Cell ', 'other_assay').status).to eq(422)
      expect(error_message).to match(/already exists/)
      other = add_modality('bulk', 'other_assay')
      expect(other.status).to eq(200), other.body
      id = json_body[:modality][:id]
      expect(update_modality(id, name: 'SINGLE CELL').status).to eq(422)
      expect(update_modality(id, name: 'Bulk').status).to eq(200)
      expect(modalities.map { |modality| modality[:name] }).to eq(['single cell', 'Bulk'])
    end

    it 'orders modalities by position then id' do
      add_modality('A', 'assay')
      add_modality('B', 'assay')
      add_modality('C', 'assay', position: 0)
      expect(modalities.map { |modality| modality[:name] }).to eq(%w[C A B])
    end

    it 'returns 404 when updating or deleting a missing modality' do
      expect(update_modality(9999, name: 'x').status).to eq(404)
      expect(error_message).to eq('No such modality')
      expect(delete_modality(9999).status).to eq(404)
    end

    it 'refuses to delete a modality that fields use and says how many' do
      add_modality('single cell', 'assay')
      modality_id = json_body[:modality][:id]
      add_field('Pools', ['subject', '::count'], modality_id: modality_id)
      pools = json_body[:field][:id]
      add_field('Kind', ['subject', '::distinct', 'species'], modality_id: modality_id)

      expect(delete_modality(modality_id).status).to eq(422)
      expect(error_message).to match(/used by 2/)
      expect(modalities.length).to eq(1)

      expect(delete_field(pools).status).to eq(200)
      expect(delete_modality(modality_id).status).to eq(422)
      expect(error_message).to match(/used by 1/)

      delete_field(fields.first[:id])
      expect(delete_modality(modality_id).status).to eq(200)
      expect(modalities).to eq([])
    end
  end

  describe CrossProjectSearchController do
    before do
      add_model('cps_example', 'patient', template_model_name: 'subject')
      add_model('cps_example', 'cohort', template_model_name: 'project')
      create_attribute('cps_example', 'patient', 'species')
      create_attribute('cps_example', 'cohort', 'project_status')
      create_attribute('cps_example', 'cohort', 'completed', type: 'boolean')
      add_table('cps_example', 'patient', { species: :string }, [
        { species: 'human' }, { species: 'mouse' }, { species: 'human' }
      ])
      add_table('cps_example', 'cohort', { project_status: :string, completed: :boolean }, [
        { project_status: 'Active', completed: false }
      ])
      add_model('cps_bare', 'unrelated')
      load_projects('cps_example', 'cps_bare')
    end

    it 'returns a count, a first value, a distinct list and a boolean, one result per query per project' do
      queries = [
        ['subject', '::count'],
        ['project', '::first', 'project_status'],
        ['subject', '::distinct', 'species'],
        ['project', '::first', 'completed']
      ]
      projects = results(queries: queries)

      expect(projects[:cps_example]).to eq([3, 'Active', ['human', 'mouse'], false])
      expect(projects).not_to have_key(:coprojects_template)
    end

    it 'asks each question as the requesting user, with restriction and the query timeout' do
      expect(Magma::Question).to receive(:new).
        with('cps_example', ['patient', '::count'], hash_including(restrict: true, user: an_instance_of(Etna::User), timeout: Magma.instance.config(:query_timeout))).
        and_call_original
      allow(Magma::Question).to receive(:new).with(satisfy { |name| name != 'cps_example' }, anything, anything).and_call_original

      expect(results(queries: [['subject', '::count']])[:cps_example]).to eq([3])
    end

    it 'returns null for a project that does not map the model, not an error' do
      projects = results(queries: [['subject', '::count'], ['project', '::first', 'project_status']])
      expect(projects[:cps_bare]).to eq([nil, nil])
    end

    it 'returns null when the project maps the model but lacks the attribute' do
      add_model('cps_thin', 'cohort', template_model_name: 'project')
      create_attribute('cps_thin', 'cohort', 'name')
      add_table('cps_thin', 'cohort', { name: :string }, [{ name: 'thin' }])
      load_projects('cps_example', 'cps_bare', 'cps_thin')
      projects = results(queries: [['project', '::first', 'project_status'], ['subject', '::count']])
      expect(projects[:cps_thin]).to eq([nil, nil])
    end

    it 'returns null for a project that is not loaded' do
      add_model('cps_unloaded', 'patient', template_model_name: 'subject')
      expect(results(queries: [['subject', '::count']])[:cps_unloaded]).to eq([nil])
    end

    it 'runs the queries of the stored rows' do
      add_field('Subjects', ['subject', '::count'])
      add_field('Species', ['subject', '::distinct', 'species'])
      queries = fields.map { |field| field[:query] }
      expect(results(queries: queries)[:cps_example]).to eq([3, ['human', 'mouse']])
    end

    it 'does not let one failing project break the others' do
      add_model('cps_broken', 'patient', template_model_name: 'subject')
      create_attribute('cps_broken', 'patient', 'species')
      add_table('cps_broken', 'patient', { species: :string }, [{ species: 'human' }])
      load_projects('cps_example', 'cps_bare', 'cps_broken')
      Magma.instance.db.drop_table?(Sequel[:cps_broken][:patients])
      # Examples run inside one transaction; with a timeout each question runs
      # in its own nested transaction, so the failed one does not abort the rest.
      allow(Magma.instance).to receive(:config).and_call_original
      allow(Magma.instance).to receive(:config).with(:query_timeout).and_return(10_000)

      expect(Magma.instance.logger).to receive(:log_error).with(an_instance_of(Sequel::DatabaseError))
      projects = results(queries: [['subject', '::count']])
      expect(projects[:cps_broken]).to eq([nil])
      expect(projects[:cps_example]).to eq([3])
    end

    it 'still counts through an intermediate model between template models' do
      add_search_lineage('cps_inter', intermediate: 'timepoint')
      add_search_lineage('cps_direct')
      load_projects('cps_example', 'cps_bare', 'cps_inter', 'cps_direct')
      query = ['subject', ['biospecimen', ['sc_seq', ['::has', 'sc_seq_pool'], '::any'], '::any'], '::count']

      projects = results(queries: [query])
      expect(projects[:cps_inter]).to eq([1])
      expect(projects[:cps_direct]).to eq([1])
    end

    it 'accepts exactly the batch cap and rejects one more' do
      expect(run_queries(queries: Array.new(50) { ['subject', '::count'] }).status).to eq(200)
      expect(run_queries(queries: Array.new(51) { ['subject', '::count'] }).status).to eq(422)
    end

    it 'names the query that is not an array and rejects the batch' do
      expect(run_queries(queries: [['subject', '::count'], 'junk']).status).to eq(422)
      expect(error_message).to eq('Query 2: Query must be a non-empty JSON array whose first item is a model name, e.g. ["model_name", "::count"]')
      expect(run_queries(queries: [[]]).status).to eq(422)
      expect(error_message).to start_with('Query 1:')
    end

    it 'returns null for a query Magma cannot run' do
      [
        ['subject', '::bogus'],
        ['subject', '::count', 'species'],
        ['not_in_template', '::count'],
        ['::count'],
        ['subject', { a: 1 }]
      ].each do |query|
        expect(results(queries: [query])[:cps_example]).to eq([nil]), query.inspect
      end
    end

    it 'rejects a batch that is not a list of queries' do
      [
        {},
        { queries: [] },
        { queries: 'subject' },
        { queries: { a: 1 } },
        { query: ['subject', '::count'] }
      ].each do |body|
        response = run_queries(body)
        expect(response.status).to eq(422), "#{body} -> #{response.status} #{response.body}"
      end
    end

    it 'lets Magma judge a deeply nested query' do
      deep = ['subject', '::count']
      12.times { deep = [deep] }
      deep = ['subject', deep]
      expect(results(queries: [deep])[:cps_example]).to eq([nil])
    end
  end
end
