require 'webmock/rspec'

describe Etna::Clients::Magma::ModelsOdmXml::Exporter do
  before(:each) do
  end

  def create_model(model_name, identifier, parent = nil, parent_link_type = nil)
    models.build_model(model_name).build_template.tap do |model_template|
      model_template.identifier = identifier
      model_template.name = model_name
      model_template.parent = parent
      model_template.build_attributes.build_attribute(identifier).tap do |attr|
        attr.attribute_name = attr.name = identifier
        attr.attribute_type = Etna::Clients::Magma::AttributeType::IDENTIFIER
        attr.display_name = attr.name
      end

      unless parent_link_type.nil?
        models.build_model(parent).build_template.build_attributes.build_attribute(model_name).tap do |attr|
          attr.attribute_name = attr.name = model_name
          attr.display_name = attr.name
          attr.attribute_type = parent_link_type
          attr.link_model_name = model_name
        end

        model_template.build_attributes.build_attribute(parent).tap do |attr|
          attr.attribute_name = attr.name = parent
          attr.display_name = attr.name
          attr.attribute_type = Etna::Clients::Magma::AttributeType::PARENT
          attr.link_model_name = parent
        end
      end
    end
  end

  let(:models) { Etna::Clients::Magma::Models.new }
  let(:root_model_template) do
    create_model('project', 'name')
  end

  let(:table_child_template) do
    create_model('table_child', 'some-id', root_model_template.name, Etna::Clients::Magma::AttributeType::TABLE)
  end

  let(:non_table_child_template) do
    create_model('non_table_child', 'some_id', root_model_template.name, Etna::Clients::Magma::AttributeType::CHILD)
  end

  it 'writes an XML file' do
    subject = create_model('subject', 'name', root_model_template.name, Etna::Clients::Magma::AttributeType::COLLECTION)
    timepoint = create_model('timepoint', 'name', subject.name, Etna::Clients::Magma::AttributeType::COLLECTION)

    project = Etna::Clients::Magma::ModelsOdmXml::Exporter.new(project_name: PROJECT)

    xml = project.write_models(models)
    puts xml
    expect(xml.include?("SourceSystem=\"Magma\"")).to eq(true)
    expect(xml.include?("Study OID=\"Project.Test\"")).to eq(true)
  end
end
