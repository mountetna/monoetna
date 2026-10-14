describe TemplateAuditController do
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

  def add_child_attribute(model_name, attribute_name)
    create_attribute('audit_project', model_name, attribute_name)
  end

  def create_model(project_name, model_name, template_project_name: nil, template_model_name: nil)
    Magma.instance.db[:models].insert(
      project_name: project_name,
      model_name: model_name,
      template_project_name: template_project_name,
      template_model_name: template_model_name
    )
  end

  def create_attribute(project_name, model_name, attribute_name, required_on_child_models: false, validation: nil)
    Magma.instance.db[:attributes].insert(
      project_name: project_name,
      model_name: model_name,
      attribute_name: attribute_name,
      column_name: attribute_name,
      type: 'string',
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

  def audit_report
    auth_header(:superuser)
    get('/template_audit')

    expect(last_response.status).to eq(200)
    json_body[:projects].find { |project| project[:project] == 'audit_project' }
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

    expect(audit_report).to eq(
      project: 'audit_project',
      conforming: true,
      unmapped_models: [],
      invalid_mappings: [],
      missing_template_columns: [],
      invalid_ontology_values: []
    )
  end

  it 'reports unmapped models and invalid mappings' do
    add_child_model('assay', template_model_name: 'unknown_model')
    add_child_model('patient', template_project_name: 'another_template', template_model_name: 'sample')
    add_child_model('subject')
    add_child_model('visit', template_project_name: 'coprojects_template')

    report = audit_report

    expect(report[:conforming]).to eq(false)
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
    expect(report[:missing_template_columns]).to eq([
      { model: 'local_sample', template_model: 'sample', columns: ['species'] }
    ])
  end

  it 'reports stored values that are not in the ontology table' do
    add_template_attribute('sample', 'tissue', validation: uberon_validation)
    uberon = double('uberon', attributes: { name: double(column_name: 'name') }, select_map: ['blood', 'lung'])
    allow(Magma.instance).to receive(:get_model).and_call_original
    allow(Magma.instance).to receive(:get_model).with('ontologies', 'uberon').and_return(uberon)

    add_child_model('local_sample', template_model_name: 'sample')
    add_child_attribute('local_sample', 'species')
    add_child_attribute('local_sample', 'tissue')
    add_child_records('local_sample', :tissue, ['blood', 'lung', 'PBMC', nil])

    report = audit_report

    expect(report[:conforming]).to eq(false)
    expect(report[:invalid_ontology_values]).to eq([
      {
        model: 'local_sample',
        template_model: 'sample',
        column: 'tissue',
        ontology: { project_name: 'ontologies', model_name: 'uberon', attribute_name: 'name' },
        values: ['PBMC']
      }
    ])
  end

  it 'excludes the template and ontology projects' do
    add_template_attribute('sample', 'tissue', validation: uberon_validation)
    create_model('ontologies', 'uberon')

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
