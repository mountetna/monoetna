describe CrossProjectSearchAuditController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  def add_template_model(model_name)
    create_model('coprojects_template', model_name)
  end

  def add_template_attribute(model_name, attribute_name, **options)
    create_attribute('coprojects_template', model_name, attribute_name, **options)
  end

  def add_child_model(model_name, template_model_name: nil, template_project_name: template_model_name && 'coprojects_template')
    create_model('audit_project', model_name, template_project_name: template_project_name, template_model_name: template_model_name)
  end

  def add_child_attribute(model_name, attribute_name, **options)
    create_attribute('audit_project', model_name, attribute_name, **options)
  end

  def create_model(project_name, model_name, template_project_name: nil, template_model_name: nil)
    Magma.instance.db[:models].insert(
      project_name: project_name,
      model_name: model_name,
      template_project_name: template_project_name,
      template_model_name: template_model_name
    )
  end

  def create_attribute(project_name, model_name, attribute_name, required_on_child_models: false, validation: nil, type: 'string', link_model_name: nil, column_name: attribute_name)
    Magma.instance.db[:attributes].insert(
      project_name: project_name,
      model_name: model_name,
      attribute_name: attribute_name,
      column_name: column_name,
      type: type,
      link_model_name: link_model_name,
      required_on_child_models: required_on_child_models,
      validation: validation && Sequel.pg_json_wrap(validation)
    )
  end

  def add_child_records(model_name, column, values)
    table = Sequel[:audit_project][model_name.pluralize.to_sym]
    Magma.instance.db.create_schema(:audit_project, if_not_exists: true)
    Magma.instance.db.create_table(table) { String column }
    values.each { |value| Magma.instance.db[table].insert(column => value) }
  end

  def add_child_foreign_key(model_name, attribute_name, type, link_model_name)
    add_child_attribute(model_name, attribute_name, type: type, link_model_name: link_model_name, column_name: "#{attribute_name}_id")
  end

  def add_child_table(model_name, columns, rows = [])
    table = Sequel[:audit_project][model_name.pluralize.to_sym]
    Magma.instance.db.create_schema(:audit_project, if_not_exists: true)
    Magma.instance.db.create_table(table) do
      primary_key :id
      columns.each { |column| Integer column }
    end
    rows.each { |row| Magma.instance.db[table].insert(columns.zip(row).to_h) }
  end

  def project_audit
    auth_header(:superuser)
    get('/cross_project_search_audit/audit_project')

    expect(last_response.status).to eq(200)
    json_body
  end

  def audit_report
    project_audit[:template_mapping]
  end

  def single_cell_report
    project_audit[:model_relationships][:single_cell]
  end

  let(:uberon_validation) do
    { type: 'Ontology', project_name: 'ontologies', model_name: 'uberon', attribute_name: 'name' }
  end

  before do
    add_template_model('sample')
    add_template_attribute('sample', 'species', required_on_child_models: true)
  end

  it 'reports a conforming project' do
    add_child_model('local_sample', template_model_name: 'sample')
    add_child_attribute('local_sample', 'species')

    report = project_audit

    expect(report).to include(project: 'audit_project', conforming: true)
    expect(report[:template_mapping]).to eq(
      conforming: true,
      summary: {
        models: { total: 1, mapped: 1, unmapped: 0, invalid: 0 },
        required_columns: { total: 1, confirmed: 1, missing: 0 },
        ontologies: { total: 0, confirmed: 0, missing_columns: 0, mismatched: 0, invalid_values: 0 }
      },
      mapped_models: [
        { model: 'local_sample', template_model: 'sample' }
      ],
      unmapped_models: [],
      invalid_mappings: [],
      missing_template_columns: [],
      ontology_mismatches: [],
      invalid_ontology_values: [],
      required_column_results: [
        { model: 'local_sample', template_model: 'sample', column: 'species', present: true }
      ],
      ontology_value_results: [],
      missing_ontology_columns: []
    )
  end

  it 'reports unmapped models and invalid mappings' do
    add_child_model('assay', template_model_name: 'unknown_model')
    add_child_model('patient', template_project_name: 'another_template', template_model_name: 'sample')
    add_child_model('project')
    add_child_model('subject')
    add_child_model('visit', template_project_name: 'coprojects_template')

    report = audit_report

    expect(report[:conforming]).to eq(false)
    expect(report[:summary][:models]).to eq(total: 4, mapped: 0, unmapped: 1, invalid: 3)
    expect(report[:mapped_models]).to eq([])
    expect(report[:unmapped_models]).to eq(['subject'])
    expect(report[:invalid_mappings]).to eq([
      { model: 'assay', template_model: 'unknown_model' },
      { model: 'patient', template_model: 'sample' },
      { model: 'visit', template_model: nil }
    ])
  end

  it 'reports missing required template columns' do
    add_child_model('local_sample', template_model_name: 'sample')

    report = audit_report

    expect(report[:conforming]).to eq(false)
    expect(report[:summary][:required_columns]).to eq(total: 1, confirmed: 0, missing: 1)
    expect(report[:missing_template_columns]).to eq([
      { model: 'local_sample', template_model: 'sample', columns: ['species'] }
    ])
    expect(report[:required_column_results]).to eq([
      { model: 'local_sample', template_model: 'sample', column: 'species', present: false }
    ])
  end

  it 'reports required columns with a missing or different ontology' do
    add_template_attribute('sample', 'tissue', required_on_child_models: true, validation: uberon_validation)

    add_child_model('blood_sample', template_model_name: 'sample')
    add_child_attribute('blood_sample', 'species')
    add_child_attribute('blood_sample', 'tissue')
    add_child_records('blood_sample', :tissue, [])

    add_child_model('local_sample', template_model_name: 'sample')
    add_child_attribute('local_sample', 'species')
    add_child_attribute('local_sample', 'tissue', validation: uberon_validation.merge(model_name: 'my_tissues'))
    add_child_records('local_sample', :tissue, [])

    report = audit_report

    expect(report[:conforming]).to eq(false)
    expect(report[:summary][:ontologies]).to eq(
      total: 2,
      confirmed: 0,
      missing_columns: 0,
      mismatched: 2,
      invalid_values: 0
    )
    expect(report[:ontology_mismatches]).to eq([
      {
        model: 'blood_sample',
        template_model: 'sample',
        column: 'tissue',
        expected: { project_name: 'ontologies', model_name: 'uberon', attribute_name: 'name' },
        actual: nil
      },
      {
        model: 'local_sample',
        template_model: 'sample',
        column: 'tissue',
        expected: { project_name: 'ontologies', model_name: 'uberon', attribute_name: 'name' },
        actual: { project_name: 'ontologies', model_name: 'my_tissues', attribute_name: 'name' }
      }
    ])
  end

  it 'reports stored values that are not in the ontology table' do
    add_template_attribute('sample', 'tissue', required_on_child_models: true, validation: uberon_validation)
    uberon = double('uberon', attributes: { name: double(column_name: 'name') }, select_map: ['blood', 'lung'])
    allow(Magma.instance).to receive(:get_model).and_call_original
    allow(Magma.instance).to receive(:get_model).with('ontologies', 'uberon').and_return(uberon)

    add_child_model('local_sample', template_model_name: 'sample')
    add_child_attribute('local_sample', 'species')
    add_child_attribute('local_sample', 'tissue', validation: uberon_validation)
    add_child_records('local_sample', :tissue, ['blood', 'lung', 'PBMC', nil])

    report = audit_report

    expect(report[:conforming]).to eq(false)
    expect(report[:summary][:ontologies]).to eq(
      total: 1,
      confirmed: 1,
      missing_columns: 0,
      mismatched: 0,
      invalid_values: 1
    )
    expect(report[:ontology_mismatches]).to eq([])
    expect(report[:invalid_ontology_values]).to eq([
      {
        model: 'local_sample',
        template_model: 'sample',
        column: 'tissue',
        ontology: { project_name: 'ontologies', model_name: 'uberon', attribute_name: 'name' },
        values: ['PBMC']
      }
    ])
    expect(report[:required_column_results]).to include(
      { model: 'local_sample', template_model: 'sample', column: 'species', present: true },
      { model: 'local_sample', template_model: 'sample', column: 'tissue', present: true }
    )
    expect(report[:ontology_value_results]).to include(
      { model: 'local_sample', template_model: 'sample', column: 'tissue', value: 'blood', valid: true },
      { model: 'local_sample', template_model: 'sample', column: 'tissue', value: 'lung', valid: true },
      { model: 'local_sample', template_model: 'sample', column: 'tissue', value: 'PBMC', valid: false }
    )
    expect(report[:missing_ontology_columns]).to eq([])
  end

  it 'lists a missing required ontology column apart from validation mismatches' do
    add_template_attribute('sample', 'tissue', required_on_child_models: true, validation: uberon_validation)

    add_child_model('local_sample', template_model_name: 'sample')
    add_child_attribute('local_sample', 'species')

    report = audit_report

    expect(report[:required_column_results]).to eq([
      { model: 'local_sample', template_model: 'sample', column: 'species', present: true },
      { model: 'local_sample', template_model: 'sample', column: 'tissue', present: false }
    ])
    expect(report[:missing_ontology_columns]).to eq([
      { model: 'local_sample', template_model: 'sample', column: 'tissue' }
    ])
    expect(report[:ontology_mismatches]).to eq([])
    expect(report[:ontology_value_results]).to eq([])
  end

  it 'ignores ontologies on columns that are not required' do
    add_template_attribute('sample', 'tissue', validation: uberon_validation)

    add_child_model('local_sample', template_model_name: 'sample')
    add_child_attribute('local_sample', 'species')
    add_child_attribute('local_sample', 'tissue')

    expect(audit_report[:conforming]).to eq(true)
  end

  it 'audits only the requested project' do
    add_child_model('local_sample', template_model_name: 'sample')
    add_child_attribute('local_sample', 'species')
    create_model('other_project', 'other_sample')

    report = project_audit

    expect(report[:project]).to eq('audit_project')
    expect(report[:template_mapping][:summary][:models]).to eq(total: 1, mapped: 1, unmapped: 0, invalid: 0)
  end

  it 'allows a project administrator to audit their project' do
    create_model('labors', 'local_sample', template_project_name: 'coprojects_template', template_model_name: 'sample')

    auth_header(:admin)
    get('/cross_project_search_audit/labors')

    expect(last_response.status).to eq(200)
    expect(json_body[:project]).to eq('labors')
  end

  it 'reports a requested project that has no models to audit' do
    auth_header(:superuser)
    get('/cross_project_search_audit/empty_project')

    expect(last_response.status).to eq(200)
    expect(json_body).to include(project: 'empty_project', conforming: true)
    expect(json_body[:template_mapping]).to include(
      conforming: true,
      summary: {
        models: { total: 0, mapped: 0, unmapped: 0, invalid: 0 },
        required_columns: { total: 0, confirmed: 0, missing: 0 },
        ontologies: { total: 0, confirmed: 0, missing_columns: 0, mismatched: 0, invalid_values: 0 }
      }
    )
  end

  it 'does not allow a project viewer to audit their project' do
    auth_header(:viewer)
    get('/cross_project_search_audit/labors')

    expect(last_response.status).to eq(403)
  end

  context 'model hierarchy' do
    before do
      add_template_model('subject')
      add_template_model('biospecimen')
      add_template_attribute('biospecimen', 'subject', type: 'parent', link_model_name: 'subject', column_name: 'subject_id')
      add_child_model('patient', template_model_name: 'subject')
      add_child_model('specimen', template_model_name: 'biospecimen')
    end

    it 'conforms when child parents follow the template' do
      add_child_foreign_key('specimen', 'patient', 'parent', 'patient')

      expect(project_audit[:model_relationships][:hierarchy]).to eq(conforming: true, issues: [])
    end

    it 'reports a parent that does not follow the template' do
      add_template_model('timepoint')
      add_child_model('visit', template_model_name: 'timepoint')
      add_child_foreign_key('specimen', 'visit', 'parent', 'visit')

      report = project_audit

      expect(report[:conforming]).to eq(false)
      expect(report[:template_mapping][:conforming]).to eq(true)
      expect(report[:model_relationships][:hierarchy]).to eq(
        conforming: false,
        issues: [
          { model: 'specimen', template_model: 'biospecimen', parent: 'visit', expected_parent: 'patient' }
        ]
      )
    end
  end

  context 'single-cell depooling' do
    before do
      %w(subject biospecimen sc_seq sc_seq_pool).each { |model_name| add_template_model(model_name) }
      add_template_attribute('biospecimen', 'subject', type: 'parent', link_model_name: 'subject', column_name: 'subject_id')
      add_template_attribute('sc_seq', 'biospecimen', type: 'parent', link_model_name: 'biospecimen', column_name: 'biospecimen_id')

      add_child_model('patient', template_model_name: 'subject')
      add_child_model('specimen', template_model_name: 'biospecimen')
      add_child_model('tube', template_model_name: 'sc_seq')
      add_child_foreign_key('specimen', 'patient', 'parent', 'patient')
      add_child_foreign_key('tube', 'specimen', 'parent', 'specimen')
    end

    it 'is not applicable when sc_seq_pool is not mapped' do
      add_child_model('pool')
      add_child_table('patient', [])
      add_child_table('specimen', [:patient_id])
      add_child_table('tube', [:specimen_id])
      add_child_table('pool', [])

      report = project_audit

      expect(report[:model_relationships][:single_cell]).to eq(
        applicable: false,
        conforming: true,
        reason: 'Child project has no model for sc_seq_pool',
        child_models: { subject: 'patient', biospecimen: 'specimen', sc_seq: 'tube', sc_seq_pool: nil }
      )
      expect(report[:conforming]).to eq(false)
      expect(report[:template_mapping][:unmapped_models]).to eq(['pool'])
    end

    it 'reports a missing link from sc_seq to sc_seq_pool' do
      add_child_model('pool', template_model_name: 'sc_seq_pool')
      add_child_table('specimen', [:patient_id])
      add_child_table('tube', [:specimen_id], [[1], [nil]])
      add_child_table('pool', [], [[], []])

      report = project_audit
      single_cell = report[:model_relationships][:single_cell]

      expect(report[:conforming]).to eq(false)
      expect(report[:template_mapping][:conforming]).to eq(true)
      expect(single_cell).to include(
        applicable: true,
        conforming: false,
        pool_link: nil,
        subject_path: ['tube', 'specimen', 'patient'],
        issues: ['tube has no link attribute to pool']
      )
      expect(single_cell[:counts]).to eq(
        pools: { total: 2, depooled: nil },
        sc_seq: { total: 2, pool_linked: nil, depooled: nil },
        subjects: nil
      )
    end

    it 'counts depooled pools, sc_seq records and subjects' do
      add_child_model('pool', template_model_name: 'sc_seq_pool')
      add_child_foreign_key('tube', 'pool', 'link', 'pool')

      add_child_table('patient', [], [[], [], []])
      add_child_table('specimen', [:patient_id], [[1], [2], [nil]])
      add_child_table('pool', [], [[], [], []])
      add_child_table('tube', [:specimen_id, :pool_id], [
        [1, 1],
        [2, 1],
        [3, 2],
        [nil, 2],
        [1, nil]
      ])

      report = project_audit

      expect(report[:conforming]).to eq(true)
      expect(report[:model_relationships][:single_cell]).to eq(
        applicable: true,
        conforming: true,
        child_models: { subject: 'patient', biospecimen: 'specimen', sc_seq: 'tube', sc_seq_pool: 'pool' },
        pool_link: 'pool',
        subject_path: ['tube', 'specimen', 'patient'],
        issues: [],
        counts: {
          pools: { total: 3, depooled: 2 },
          sc_seq: { total: 5, pool_linked: 4, depooled: 3 },
          subjects: 2
        }
      )
    end
  end
end
