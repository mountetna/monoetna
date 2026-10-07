describe TemplateAuditController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  def add_model(project_name, model_name, template_model_name: nil, template_project_name: template_model_name && 'coprojects_template')
    Magma.instance.db[:models].insert(
      project_name: project_name,
      model_name: model_name,
      template_project_name: template_project_name,
      template_model_name: template_model_name
    )
  end

  def add_attribute(project_name, model_name, attribute_name, template_required: false, validation: nil)
    Magma.instance.db[:attributes].insert(
      project_name: project_name,
      model_name: model_name,
      attribute_name: attribute_name,
      column_name: attribute_name,
      type: 'string',
      template_required: template_required,
      validation: validation && Sequel.pg_json_wrap(validation)
    )
  end

  def add_records(project_name, model_name, column, values)
    table = Sequel[project_name.to_sym][model_name.pluralize.to_sym]
    Magma.instance.db.create_schema(project_name.to_sym, if_not_exists: true)
    Magma.instance.db.create_table(table) { String column }
    values.each { |value| Magma.instance.db[table].insert(column => value) }
  end

  def audit_report(project_name)
    auth_header(:superuser)
    get('/template_audit')

    expect(last_response.status).to eq(200)
    json_body[:projects].find { |project| project[:project] == project_name }
  end

  before do
    add_model('coprojects_template', 'sample')
    add_attribute('coprojects_template', 'sample', 'species', template_required: true)
  end

  it 'reports a conforming project' do
    add_model('audit_project', 'local_sample', template_model_name: 'sample')
    add_attribute('audit_project', 'local_sample', 'species')

    expect(audit_report('audit_project')).to eq(
      project: 'audit_project',
      conforming: true,
      unmapped_models: [],
      invalid_mappings: [],
      missing_template_columns: [],
      invalid_ontology_values: []
    )
  end

  it 'reports unmapped models and invalid mappings' do
    add_model('audit_project', 'assay', template_model_name: 'unknown_model')
    add_model('audit_project', 'patient', template_project_name: 'another_template', template_model_name: 'sample')
    add_model('audit_project', 'subject')
    add_model('audit_project', 'visit', template_project_name: 'coprojects_template')

    report = audit_report('audit_project')

    expect(report[:conforming]).to eq(false)
    expect(report[:unmapped_models]).to eq(['subject'])
    expect(report[:invalid_mappings]).to eq([
      { model: 'assay', template_model: 'unknown_model' },
      { model: 'patient', template_model: 'sample' },
      { model: 'visit', template_model: nil }
    ])
  end

  it 'reports missing required template columns' do
    add_model('audit_project', 'local_sample', template_model_name: 'sample')

    report = audit_report('audit_project')

    expect(report[:conforming]).to eq(false)
    expect(report[:missing_template_columns]).to eq([
      { model: 'local_sample', template_model: 'sample', columns: ['species'] }
    ])
  end

  it 'reports stored values that are not in the ontology table' do
    add_attribute('coprojects_template', 'sample', 'tissue', validation: { type: 'Ontology', value: 'uberon' })
    uberon = double('uberon', identity: double(column_name: :name), select_map: ['blood', 'lung'])
    allow(Magma.instance).to receive(:get_model).and_call_original
    allow(Magma.instance).to receive(:get_model).with('ontologies', 'uberon').and_return(uberon)

    add_model('audit_project', 'local_sample', template_model_name: 'sample')
    add_attribute('audit_project', 'local_sample', 'species')
    add_attribute('audit_project', 'local_sample', 'tissue')
    add_records('audit_project', 'local_sample', :tissue, ['blood', 'lung', 'PBMC', nil])

    report = audit_report('audit_project')

    expect(report[:conforming]).to eq(false)
    expect(report[:invalid_ontology_values]).to eq([
      { model: 'local_sample', template_model: 'sample', column: 'tissue', table: 'uberon', values: ['PBMC'] }
    ])
  end

  it 'excludes the template and ontologies projects' do
    add_model('ontologies', 'uberon')

    auth_header(:superuser)
    get('/template_audit')

    expect(json_body[:projects].map { |project| project[:project] }).not_to include('coprojects_template', 'ontologies')
  end

  it 'requires a supereditor' do
    auth_header(:viewer)
    get('/template_audit')

    expect(last_response.status).to eq(403)
  end
end
