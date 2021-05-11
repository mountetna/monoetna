require 'nokogiri'

module Etna
  module Clients
    class Magma
      module ModelsOdmXml
      end
    end
  end
end

module Etna
  module Clients
    class Magma
      module ModelsOdmXml
        module Prettify
          def prettify(name)
            name.split('_').map(&:capitalize).join(' ')
          end

          def shorten(name)
            name.gsub('_', '').capitalize
          end
        end

        class Exporter
          include Prettify

          attr_reader :project_name, :models

          def initialize(project_name:, models:)
            @project_name = project_name
            @models = models
          end

          def data_type_map
            @data_type_map ||= begin
              map = {}
              map[Etna::Clients::Magma::AttributeType::STRING] = 'text'
              map[Etna::Clients::Magma::AttributeType::DATE_TIME] = 'date'
              map[Etna::Clients::Magma::AttributeType::BOOLEAN] = 'text'
              map[Etna::Clients::Magma::AttributeType::FLOAT] = 'float'
              map[Etna::Clients::Magma::AttributeType::INTEGER] = 'integer'

              map
            end
          end

          def redcap_field_type_map
            @redcap_field_type_map ||= begin
              map = {}
              map[Etna::Clients::Magma::AttributeType::STRING] = 'textarea'
              map[Etna::Clients::Magma::AttributeType::DATE_TIME] = 'text'
              map[Etna::Clients::Magma::AttributeType::BOOLEAN] = 'radio'
              map[Etna::Clients::Magma::AttributeType::FLOAT] = 'text'
              map[Etna::Clients::Magma::AttributeType::INTEGER] = 'text'

              map
            end
          end

          def redcap_text_validation_map
            @redcap_text_validation_map ||= begin
              map = {}
              map[Etna::Clients::Magma::AttributeType::DATE_TIME] = 'date_mdy'
              map[Etna::Clients::Magma::AttributeType::FLOAT] = 'float'
              map[Etna::Clients::Magma::AttributeType::INTEGER] = 'int'

              map
            end
          end

          def write_models(output_io: nil, filename: nil)
            @document = Nokogiri::XML::Builder.new(encoding: 'UTF-8') do |xml|
              xml.ODM(odm_headers) do
                xml.Study(OID: "Project.#{shorten(project_name)}") do
                end
                global_variables(xml)
                metadata(xml)
              end
            end

            @document.to_xml
          end

          def odm_headers
            {
              xmlns: "http://www.cdisc.org/ns/odm/v1.3",
              'xmlns:ds': "http://www.w3.org/2000/09/xmldsig#",
              'xmlns:xsi': "http://www.w3.org/2001/XMLSchema-instance",
              'xmlns:redcap': "https://projectredcap.org",
              'xsi:schemaLocation': "http://www.cdisc.org/ns/odm/v1.3 schema/odm/ODM1-3-1.xsd",
              ODMVersion: "1.3.2",
              FileOID: "000-00-0000",
              FileType: "Snapshot",
              Description: project_name,
              AsOfDateTime: DateTime.now,
              CreationDateTime: DateTime.now,
              SourceSystem: "Magma",
              SourceSystemVersion: DateTime.now
            }
          end

          def global_variables(xml)
            # Includes general metadata about the project, as well as
            # declarations of all repeating instruments,
            #   which seem like Timepoints to me.
            # NOTE:
            # <redcap:Purpose>0</redcap:Purpose>
            # 0 = Practice / just for fun
            # 1 = Operational Support
            # 2 = Research
            # 3 = Quality Improvement
            # 4 = Other

            xml.GlobalVariables do
              xml.StudyName "#{project_name}"
              xml.StudyDescription "#{project_name} - Data Library integration project"
              xml.ProtocolName "#{project_name}"
              xml.send('redcap:RecordAutonumberingEnabled', 1)
              xml.send('redcap:CustomRecordLabel')
              xml.send('redcap:SecondaryUniqueField')
	            xml.send('redcap:SchedulingEnabled', 0)
	            xml.send('redcap:SurveysEnabled', 0)
              xml.send('redcap:SurveyInvitationEmailField')
              xml.send('redcap:Purpose', 2) # 2 == research
              xml.send('redcap:PurposeOther')
              xml.send('redcap:ProjectNotes', "Used to easily ingest clinical data for #{project_name} into the Data Library.")
              xml.send('redcap:MissingDataCodes')

              repeating_instruments(xml)
            end
          end

          def repeating_attribute_types
            # We don't have a good indicator for what is a repeating
            #   attribute for REDCap...
            # Will this work for models without timepoint, that go
            #   straight from subject -> sample?
            models.model('subject').template.attributes.all.select do |attribute|
              Etna::Clients::Magma::AttributeType::COLLECTION == attribute.attribute_type
            end
          end

          def clinical_dictionaries
            # Our indicator for if something needs a REDCap form will be any
            #   model with a dictionary.
            models.all.select do |model|
              model.template.raw['dictionary']
            end.map do |model|
              model.template.dictionary
            end
          end

          def repeating_instruments(xml)
            # Now we get into repeating instruments and events.
            # From a Magma model perspective, this should be
            #   Timepoint that hangs off of
            #   a Subject model.
            xml.send('redcap:RepeatingInstrumentsAndEvents') do
              repeating_attribute_types.map do |repeating_attribute|
                write_repeating_instrument_xml(xml, repeating_attribute)
              end
            end
          end

          def write_repeating_instrument_xml(xml, instrument)
            node = xml.send('redcap:RepeatingInstrument')
            node['redcap:UniqueEventName'] = 'event_1_arm_1'
            node['redcap:RepeatInstrument'] = instrument.attribute_name
            node['redcap:CustomLabel'] = instrument.display_name || instrument.attribute_name.capitalize

            node
          end

          def metadata(xml)
            # Includes form and field definitions
            xml.MetaDataVersion(
              OID: "Metadata.#{shorten(project_name)}_#{DateTime.now}",
              Name: project_name,
              'redcap:RecordIdField': 'record_id'
            ) do
              clinical_dictionaries.map do |dictionary|
                # Each Magma dictionary needs a FormDef, with
                #   ItemGroupRef children for each form page (?).
                # Each ItemGroupRef requires a corresponding ItemGroupDef
                #   with ItemRef children for each input (?).
                # Each ItemRef requires a correspdonding ItemDef
                #   that defines the label and type, and includes
                #   <Question> as a label (?).
                # Option validations are present as a CodeList (with
                #   CodeListItem children).
                write_form_def(xml, dictionary)
                write_item_group_def(xml, dictionary)
                write_item_def(xml, dictionary)
                write_code_list(xml, dictionary)
              end
            end
          end

          def write_form_def(xml, dictionary)
            xml.FormDef(
              OID: "Form.#{dictionary.model_name}",
              Name: dictionary.model_name.capitalize,
              Repeating: "No",
              'redcap:FormName': dictionary.model_name
            ) do
              # Assume a single item group
              xml.ItemGroupRef(
                ItemGroupOID: "#{dictionary.model_name}.attributes",
                Mandatory: "No"
              ) do
              end
            end
          end

          def write_item_group_def(xml, dictionary)
            xml.ItemGroupDef(
              OID: "#{dictionary.model_name}.attributes",
              Name: "#{dictionary.model_name.capitalize} Attributes",
              Repeating: "No"
            ) do
              dictionary.attributes.keys.map do |attribute_name|
                xml.ItemRef(
                  ItemOID: "#{dictionary.model_name}.#{attribute_name}", # Does this need to be unique across all items?
                  Mandatory: "No",
                  'redcap:Variable': attribute_name # Does this need to be unique?
                ) do
                end
              end
            end
          end

          def write_item_def(xml, dictionary)
            model_attributes = models.model(dictionary.model_name).template.attributes
            dictionary.attributes.keys.map do |attribute_name|
              attribute = model_attributes.attribute(attribute_name)
              attribute_type = attribute.attribute_type
              params = {
                OID: "#{dictionary.model_name}.#{attribute_name}", # Does this need to be unique across all items?
                Name: attribute_name,
                DataType: data_type_map[attribute_type] || 'text',
                'redcap:Variable': attribute_name,
                'redcap:FieldType': redcap_field_type_map[attribute_type] || 'text',
                Length: '999' # How could we infer shorter values?
              }

              params['redcap:TextValidationType'] = redcap_text_validation_map[attribute_type] if redcap_text_validation_map[attribute_type]
              params['redcap:FieldNote'] = attribute.description if attribute.description
              xml.ItemDef(params) do
                xml.Question do
                  xml.send('TranslatedText', attribute_name.capitalize)
                end

                if attribute.validation && Etna::Clients::Magma::AttributeValidationType::ARRAY == attribute.validation['type']
                  xml.CodeListRef(
                    CodeListOID: "#{dictionary.model_name}.#{attribute_name}.choices"
                  ) do
                  end
                end
              end
            end
          end

          def write_code_list(xml, dictionary)
            model_attributes = models.model(dictionary.model_name).template.attributes
            dictionary.attributes.keys.map do |attribute_name|
              attribute = model_attributes.attribute(attribute_name)
              if attribute.validation && Etna::Clients::Magma::AttributeValidationType::ARRAY == attribute.validation['type']
                xml.CodeList(
                  OID: "#{dictionary.model_name}.#{attribute_name}.choices",
                  Name: "#{attribute_name}",
                  DataType: "text",
                  'redcap:Variable': attribute_name
                ) do
                  attribute.validation['value'].map do |option|
                    xml.CodeListItem(
                      CodedValue: option
                    ) do
                      xml.Decode() do
                        xml.send('TranslatedText', option)
                      end
                    end
                  end
                end
              end
            end
          end
        end
      end
    end
  end
end