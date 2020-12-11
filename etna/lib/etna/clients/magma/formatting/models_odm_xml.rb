require 'nokogiri'

module Etna
  module Clients
    class Magma
      module ModelsOdmXml
      end
    end
  end
end

# NOTE:
# <redcap:Purpose>0</redcap:Purpose>
# 0 = Practice / just for fun
# 1 = Operational Support
# 2 = Research
# 3 = Quality Improvement
# 4 = Other

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

          attr_reader :project_name

          def initialize(project_name:)
            @project_name = project_name
          end

          def write_models(models, output_io: nil, filename: nil)
            @document = Nokogiri::XML::Builder.new(encoding: 'UTF-8') do |xml|
              xml.ODM(odm_headers) do
                xml.Study(OID: "Project.#{shorten(project_name)}") do

                end
                global_variables(xml, models)
                metadata(xml, models)
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

          def global_variables(xml, models)
            # Includes general metadata about the project, as well as
            # declarations of all repeating instruments,
            #   which seem like Timepoints to me.
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

              repeating_instruments(xml, models)
            end
          end

          def repeating_attribute_types(models)
            # Will this work for models without timepoint, that go
            #   straight from subject -> sample?
            models.model('subject').template.attributes.all.select do |attribute|
              Etna::Clients::Magma::AttributeType::COLLECTION == attribute.attribute_type
            end
          end

          def clinical_attribute_types(models)
            models.model('subject').template.attributes.all.select do |attribute|
              Etna::Clients::Magma::AttributeType::TABLE == attribute.attribute_type
            end
          end

          def repeating_instruments(xml, models)
            # Now we get into repeating instruments and events.
            # From a Magma model perspective, this should be
            #   Timepoint that hangs off of
            #   a Subject model.
            xml.send('redcap:RepeatingInstrumentsAndEvents') do
              repeating_attribute_types(models).each do |repeating_attribute|
                node = xml.send('redcap:RepeatingInstrument')
                node['redcap:UniqueEventName'] = 'event_1_arm_1'
                node['redcap:RepeatInstrument'] = repeating_attribute.attribute_name
                node['redcap:CustomLabel'] = repeating_attribute.display_name || repeating_attribute.attribute_name.capitalize

                node
              end
            end
          end

          def metadata(xml, models)
            # Includes form and field definitions
            xml.MetaDataVersion(
              OID: "Metadata.#{shorten(project_name)}_#{DateTime.now}"
              Name: project_name,
              'redcap:RecordIdField': 'record_id'
            ) do
              # Each Magma model we want here needs a FormDef, with
              #   ItemGroupRef children for each rorm page (?).
              # Each ItemGroupRef requires a corresponding ItemGroupDef
              #   with ItemRef children for each input (?).
              # Each ItemRef requires a correspdonding ItemDef
              #   that defines the label and type. and includes
              #   <Question> as a label (?).
              # Option validations are present as a CodeList (with
              #   CodeListItem children).
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

              repeating_instruments(xml, models)
            end
          end
        end
      end
    end
  end
end
