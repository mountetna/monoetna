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

  it 'generates an XML doc with ODM metadata' do
    subject = create_model('subject', 'name', root_model_template.name, Etna::Clients::Magma::AttributeType::COLLECTION)
    timepoint = create_model('timepoint', 'name', subject.name, Etna::Clients::Magma::AttributeType::COLLECTION)

    project = Etna::Clients::Magma::ModelsOdmXml::Exporter.new(
      project_name: PROJECT,
      models: models)

    xml = project.write_models

    expect(xml.include?("SourceSystem=\"Magma\"")).to eq(true)
    expect(xml.include?("Study OID=\"Project.Test\"")).to eq(true)
  end

  it 'writes forms and inputs for attributes' do
    subject = create_model('subject', 'name', root_model_template.name, Etna::Clients::Magma::AttributeType::COLLECTION)
    codex = create_model('codex', 'name', root_model_template.name, Etna::Clients::Magma::AttributeType::TABLE)

    subject = models.model('subject')
    attr_builder = subject.template.build_attributes
    date_field = attr_builder.build_attribute('favorite_date')
    date_field.attribute_type = Etna::Clients::Magma::AttributeType::DATE_TIME

    float_field = attr_builder.build_attribute('favorite_float')
    float_field.attribute_type = Etna::Clients::Magma::AttributeType::FLOAT

    integer_field = attr_builder.build_attribute('favorite_int')
    integer_field.attribute_type = Etna::Clients::Magma::AttributeType::INTEGER


    codex = models.model('codex')
    attr_builder = codex.template.build_attributes
    date_reference = attr_builder.build_attribute('date_reference')
    date_reference.attribute_type = Etna::Clients::Magma::AttributeType::DATE_TIME

    float_reference = attr_builder.build_attribute('float_reference')
    float_reference.attribute_type = Etna::Clients::Magma::AttributeType::FLOAT

    int_reference = attr_builder.build_attribute('int_reference')
    int_reference.attribute_type = Etna::Clients::Magma::AttributeType::INTEGER

    dictionary = root_model_template.build_dictionary
    dictionary.dictionary_model = "#{PROJECT}::Codex"
    dictionary.model_name = 'subject'
    dictionary.attributes = {
      'name' => 'name',
      'favorite_date' => 'date_reference',
      'favorite_float' => 'float_reference',
      'favorite_int' => 'int_reference'
    }

    project = Etna::Clients::Magma::ModelsOdmXml::Exporter.new(
      project_name: PROJECT,
      models: models)
    xml = project.write_models

    expected_items = <<-EOM
    <FormDef OID="Form.subject" Name="Subject" Repeating="No" redcap:FormName="subject">
      <ItemGroupRef ItemGroupOID="subject.attributes" Mandatory="No"/>
    </FormDef>
    <ItemGroupDef OID="subject.attributes" Name="Subject Attributes" Repeating="No">
      <ItemRef ItemOID="subject.name" Mandatory="No" redcap:Variable="name"/>
      <ItemRef ItemOID="subject.favorite_date" Mandatory="No" redcap:Variable="favorite_date"/>
      <ItemRef ItemOID="subject.favorite_float" Mandatory="No" redcap:Variable="favorite_float"/>
      <ItemRef ItemOID="subject.favorite_int" Mandatory="No" redcap:Variable="favorite_int"/>
    </ItemGroupDef>
    <ItemDef OID="subject.name" Name="name" DataType="text" redcap:Variable="name" redcap:FieldType="text" Length="999">
      <Question>
        <TranslatedText>Name</TranslatedText>
      </Question>
    </ItemDef>
    <ItemDef OID="subject.favorite_date" Name="favorite_date" DataType="date" redcap:Variable="favorite_date" redcap:FieldType="text" Length="999" redcap:TextValidationType="date_mdy">
      <Question>
        <TranslatedText>Favorite_date</TranslatedText>
      </Question>
    </ItemDef>
    <ItemDef OID="subject.favorite_float" Name="favorite_float" DataType="float" redcap:Variable="favorite_float" redcap:FieldType="text" Length="999" redcap:TextValidationType="float">
      <Question>
        <TranslatedText>Favorite_float</TranslatedText>
      </Question>
    </ItemDef>
    <ItemDef OID="subject.favorite_int" Name="favorite_int" DataType="integer" redcap:Variable="favorite_int" redcap:FieldType="text" Length="999" redcap:TextValidationType="int">
      <Question>
        <TranslatedText>Favorite_int</TranslatedText>
      </Question>
    </ItemDef>
EOM
    expect(xml.include?(expected_items)).to eq(true)
  end

  it 'writes options for attributes that have Array validation' do
    subject = create_model('subject', 'name', root_model_template.name, Etna::Clients::Magma::AttributeType::COLLECTION)
    codex = create_model('codex', 'name', root_model_template.name, Etna::Clients::Magma::AttributeType::TABLE)

    subject.attributes.attribute('name').validation = {
      'type' => 'Array',
      'value' => ['subject1', 'subject2']
    }

    dictionary = root_model_template.build_dictionary
    dictionary.dictionary_model = "#{PROJECT}::Codex"
    dictionary.model_name = 'subject'
    dictionary.attributes = {
      'name' => 'name'
    }

    project = Etna::Clients::Magma::ModelsOdmXml::Exporter.new(
      project_name: PROJECT,
      models: models)

    xml = project.write_models

    expected_code_list = <<-EOM
    <CodeList OID="subject.name.choices" Name="name" DataType="text" redcap:Variable="name">
      <CodeListItem CodedValue="subject1">
        <Decode>
          <TranslatedText>subject1</TranslatedText>
        </Decode>
      </CodeListItem>
      <CodeListItem CodedValue="subject2">
        <Decode>
          <TranslatedText>subject2</TranslatedText>
        </Decode>
      </CodeListItem>
    </CodeList>
EOM

    expect(xml.include?("<CodeListRef CodeListOID=\"subject.name.choices\"/>")).to eq(true)
    expect(xml.include?(expected_code_list)).to eq(true)
  end
end
