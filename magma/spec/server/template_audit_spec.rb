describe TemplateAuditController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  def add_model(project_name, model_name, template_project_name: nil, template_model_name: nil)
    Magma.instance.db[:models].insert(
      project_name: project_name,
      model_name: model_name,
      template_project_name: template_project_name,
      template_model_name: template_model_name
    )
  end

  def add_attribute(project_name, model_name, attribute_name, type:, template_required: false, validation: nil)
    Magma.instance.db[:attributes].insert(
      project_name: project_name,
      model_name: model_name,
      attribute_name: attribute_name,
      column_name: attribute_name,
      type: type,
      template_required: template_required,
      validation: validation
    )
  end

  def create_schema(name)
    Magma.instance.db.create_schema(name.to_sym)
  rescue Sequel::DatabaseError
    nil
  end

  def stub_ontology_table(table, terms)
    model = double(table, identity: double(column_name: :name), all: terms)

    allow(Magma.instance).to receive(:get_model).and_call_original
    allow(Magma.instance).to receive(:get_model).with('ontologies', table).and_return(model)
  end

  def set_template_validation(model_name, attribute_name, validation)
    Magma.instance.db[:attributes].
      where(
        project_name: 'coprojects_template',
        model_name: model_name,
        attribute_name: attribute_name
      ).
      update(validation: Sequel.pg_json_wrap(validation))
  end

  def create_record_table(project_name, model_name, *columns)
    create_schema(project_name)
    Magma.instance.db.create_table?(Sequel[project_name.to_sym][model_name.pluralize.to_sym]) do
      String :name
      columns.each { |column| String column.to_sym }
    end
  end

  def add_record(project_name, model_name, values)
    Magma.instance.db[Sequel[project_name.to_sym][model_name.pluralize.to_sym]].insert(values)
  end

  before do
    add_model('coprojects_template', 'sample')
    add_model('coprojects_template', 'assay')
    add_attribute(
      'coprojects_template',
      'sample',
      'species',
      type: 'string',
      template_required: true
    )
    add_attribute(
      'coprojects_template',
      'assay',
      'platform',
      type: 'string',
      template_required: true
    )
  end

  it 'reports a conforming project' do
    add_model(
      'audit_project',
      'local_sample',
      template_project_name: 'coprojects_template',
      template_model_name: 'sample'
    )
    add_model(
      'audit_project',
      'local_assay',
      template_project_name: 'coprojects_template',
      template_model_name: 'assay'
    )
    add_attribute('audit_project', 'local_sample', 'species', type: 'string')
    add_attribute('audit_project', 'local_assay', 'platform', type: 'string')

    auth_header(:superuser)
    get('/template_audit/')

    expect(last_response.status).to eq(200)
    report = json_body[:projects].find { |project| project[:project] == 'audit_project' }
    expect(report).to eq(
      project: 'audit_project',
      conforming: true,
      unmapped_models: [],
      invalid_mappings: [],
      missing_template_columns: [],
      invalid_ontology_values: []
    )
  end

  it 'reports unmapped models and invalid mappings' do
    add_model('audit_project', 'sample')
    add_model(
      'audit_project',
      'assay',
      template_project_name: 'coprojects_template',
      template_model_name: 'unknown_model'
    )

    auth_header(:superuser)
    get('/template_audit/')

    expect(last_response.status).to eq(200)
    report = json_body[:projects].find { |project| project[:project] == 'audit_project' }
    expect(report).to eq(
      project: 'audit_project',
      conforming: false,
      unmapped_models: ['sample'],
      invalid_mappings: [
        {
          model: 'assay',
          template_model: 'unknown_model'
        }
      ],
      missing_template_columns: [],
      invalid_ontology_values: []
    )
  end

  it 'treats incomplete and wrong-project references as invalid' do
    add_model(
      'audit_project',
      'assay',
      template_project_name: 'another_template',
      template_model_name: 'assay'
    )
    add_model(
      'audit_project',
      'sample',
      template_project_name: 'coprojects_template'
    )

    auth_header(:superuser)
    get('/template_audit/')

    report = json_body[:projects].find { |project| project[:project] == 'audit_project' }
    expect(report[:invalid_mappings]).to eq([
      { model: 'assay', template_model: 'assay' },
      { model: 'sample', template_model: nil }
    ])
  end

  it 'reports missing enforced template columns without comparing types' do
    add_model(
      'audit_project',
      'sample',
      template_project_name: 'coprojects_template',
      template_model_name: 'sample'
    )
    add_model(
      'audit_project',
      'assay',
      template_project_name: 'coprojects_template',
      template_model_name: 'assay'
    )
    add_attribute('audit_project', 'assay', 'platform', type: 'integer')

    auth_header(:superuser)
    get('/template_audit/')

    report = json_body[:projects].find { |project| project[:project] == 'audit_project' }
    expect(report[:conforming]).to eq(false)
    expect(report[:missing_template_columns]).to eq([
      {
        model: 'sample',
        template_model: 'sample',
        columns: ['species']
      }
    ])
    expect(report).not_to have_key(:incompatible_template_columns)
  end

  it 'reports stored values that are not in the ontology table' do
    set_template_validation('sample', 'species', type: 'Ontology', value: 'ncbitaxon')

    add_model(
      'audit_project',
      'local_sample',
      template_project_name: 'coprojects_template',
      template_model_name: 'sample'
    )
    add_attribute('audit_project', 'local_sample', 'species', type: 'string')

    stub_ontology_table('ncbitaxon', [
      { name: 'Homo sapiens', ontology_id: 'NCBITaxon:9606' },
      { name: 'Mus musculus', ontology_id: 'NCBITaxon:10090' }
    ])

    create_record_table('audit_project', 'local_sample', :species)
    add_record('audit_project', 'local_sample', name: 'S1', species: 'Homo sapiens')
    add_record('audit_project', 'local_sample', name: 'S2', species: 'NCBITaxon:10090')
    add_record('audit_project', 'local_sample', name: 'S3', species: 'Human')
    add_record('audit_project', 'local_sample', name: 'S4', species: 'Lion')
    add_record('audit_project', 'local_sample', name: 'S5', species: nil)

    auth_header(:superuser)
    get('/template_audit/')

    report = json_body[:projects].find { |project| project[:project] == 'audit_project' }
    expect(report[:conforming]).to eq(false)
    expect(report[:invalid_ontology_values]).to eq([
      {
        model: 'local_sample',
        template_model: 'sample',
        column: 'species',
        table: 'ncbitaxon',
        values: ['Human', 'Lion']
      }
    ])
  end

  it 'reports every project except the template and ontologies projects' do
    add_model('audit_project', 'sample')
    add_model('ontologies', 'uberon')

    auth_header(:superuser)
    get('/template_audit/')

    project_names = json_body[:projects].map { |project| project[:project] }
    expect(project_names).to include('audit_project', 'labors')
    expect(project_names).not_to include('coprojects_template', 'ontologies')
  end

  it 'requires a supereditor' do
    auth_header(:viewer)
    get('/template_audit/')

    expect(last_response.status).to eq(403)
  end
end
