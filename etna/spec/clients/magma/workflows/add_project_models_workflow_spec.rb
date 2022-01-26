describe Etna::Clients::Magma::AddProjectModelsWorkflow do
  describe "unit" do
    describe "downloading models into a csv" do
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
        create_model('project', 'name').tap do |root_model_template|
          root_model_template.build_attributes.build_attribute('some_num').tap do |attr|
            attr.attribute_name = attr.name = 'some_num'
            attr.display_name = attr.name
            attr.attribute_type = Etna::Clients::Magma::AttributeType::INTEGER
          end
        end
      end

      let(:table_child_template) do
        create_model('table_child', 'some-id', root_model_template.name, Etna::Clients::Magma::AttributeType::TABLE)
      end

      let(:non_table_child_template) do
        create_model('non_table_child', 'some_id', root_model_template.name, Etna::Clients::Magma::AttributeType::CHILD)
      end

      let(:root_matrix_attribute) do
        root_model_template.build_attributes.build_attribute('matrix_attr').tap do |attr|
          attr.attribute_name = attr.name = 'matrix_attr'
          attr.display_name = attr.name
          attr.attribute_type = Etna::Clients::Magma::AttributeType::MATRIX
          attr.raw['options'] = ['a', 'b', 'c']
        end
      end

      let(:root_options_attribute) do
        root_model_template.build_attributes.build_attribute('options_attr').tap do |attr|
          attr.attribute_name = attr.name = 'options_attr'
          attr.display_name = attr.name
          attr.attribute_type = Etna::Clients::Magma::AttributeType::STRING
          attr.raw['options'] = ['d', 'e', 'f']
        end
      end

      let(:rows_as_hashes) do
        result = []
        exporter = Etna::Clients::Magma::ModelsCsv::Exporter.new
        writeable = Etna::CsvExporter::RowWriteable.new(exporter, result)
        exporter.write_model_rows(models, models.model_keys.sort, writeable)
        result.map { |r| exporter.column_headers.zip(r).to_h }
      end

      describe '#validate_changeset' do
        let(:validation_errors) do
          [].tap do |errs|
            Etna::Clients::Magma::AddProjectModelsWorkflow.validate_changeset(changeset) { |err| errs << err }
          end
        end
        let(:changeset) { Etna::Clients::Magma::ModelsCsv::ModelsChangeset.new(models: models) }

        describe 'add model validations' do
          describe 'with an attribute type that does not validate' do
            before(:each) do
              bad_attr_model = create_model('bad_attr_model', 'name', root_model_template.name, Etna::Clients::Magma::AttributeType::CHILD)
              bad_attr_model.build_attributes.build_attribute('afc').tap do |attr|
                attr.name = attr.attribute_name = 'bad_attr'
                attr.attribute_type = 'not-a-thing'
              end
            end

            it 'shows that error' do
              expect(validation_errors.length).to eql(1)
              expect(validation_errors.first).to match(/Invalid attribute_type for attribute bad_attr/)
            end
          end

          describe 'with an model / attribute name that does not match validation' do
            before(:each) do
              create_model('not a good name', 'name', root_model_template.name, Etna::Clients::Magma::AttributeType::CHILD)
            end

            it 'shows that error' do
              expect(validation_errors.length).to eql(2)
              expect(validation_errors.first).to match(/Attribute name \"not a good name\"/)
              expect(validation_errors[1]).to match(/Model name not a good name must/)
            end
          end

          describe 'with identifiers / parent attributes with the wrong type' do
            before(:each) do
              create_model('with_bad_ids', 'id_is_not_id', root_model_template.name, Etna::Clients::Magma::AttributeType::CHILD).tap do |template|
                template.build_attributes.build_attribute(template.identifier).attribute_type = Etna::Clients::Magma::AttributeType::STRING
                template.build_attributes.build_attribute(template.parent).attribute_type = Etna::Clients::Magma::AttributeType::STRING
              end
            end

            it 'shows that error' do
              expect(validation_errors.length).to eql(4)
              expect(validation_errors.first).to match(/Invalid attribute_type for attribute with_bad_ids/)
              expect(validation_errors[1]).to match(/Invalid attribute_type for attribute id_is_not_id/)
              expect(validation_errors[2]).to match(/attribute project has link_model_name set, but has attribute_type/)
              expect(validation_errors[3]).to match(/Invalid attribute_type for attribute project/)
            end
          end
        end
      end

      describe '#each_csv_row' do
        describe 'table identifiers' do
          before do
            table_child_template
            non_table_child_template
          end

          it 'adds identifier rows for models, except those that are tables of parent models' do
            attribute_names = rows_as_hashes.map { |hash| hash[:attribute_name] }
            expect(attribute_names).to_not include(table_child_template.identifier)
            expect(attribute_names).to include(non_table_child_template.identifier)
          end
        end
      end

      describe 'child and parent attributes' do
        before do
          non_table_child_template
        end

        it 'is not included in the spreadsheet, as they are not configurable and are implicit' do
          attribute_names = rows_as_hashes.map { |hash| hash[:attribute_type] }
          expect(attribute_names).to_not include('parent')
          expect(attribute_names).to include('integer')
        end
      end

      describe 'matrix attributes' do
        before do
          root_matrix_attribute
          root_options_attribute
        end

        it 'includes a sentinel value instead of the actual options' do
          options = rows_as_hashes.map { |hash| hash[:options] }
          expect(options).to include('d, e, f')
          expect(options).to include(Etna::Clients::Magma::ModelsCsv::COPY_OPTIONS_SENTINEL + "64f47382e7ddc46583bf6d2abedf4140")
          expect(options).to_not include('a, b, c')
        end
      end
    end
  end

  describe "e2e" do
    it 'can sync a sub model tree from point A to point B' do
      configure_etna_yml

      VCR.use_cassette('add_project_models_workflow-timepoint-project.e2e') do
        # Change this name when re-recording the cassette file to ensure a new project is synced
        test_project = "test_add_project_models_workflow_full_ezdk"
        magma_client = Etna::Clients::Magma.new(
            host: 'https://magma.development.local',
            token: ENV['TOKEN'] || TEST_TOKEN, ignore_ssl: true,
        )

        magma_client.update_model(Etna::Clients::Magma::UpdateModelRequest.new(
            project_name: test_project,
            actions: [Etna::Clients::Magma::AddProjectAction.new]))

        # Ensure this is a brand new project we are testing with.
        expect(magma_client.retrieve(Etna::Clients::Magma::RetrievalRequest.new(project_name: test_project)).models.model_keys).to eql(['project'])

        workflow = Etna::Clients::Magma::AddProjectModelsWorkflow.new(magma_client: magma_client)

        io = StringIO.new
        workflow.write_models_template_csv('mvir1', 'timepoint', io: io)

        io.rewind

        changeset = workflow.prepare_changeset_from_csv(io: io) do |err|
          raise err
        end

        # Change this based on source model
        expect(changeset.matrix_constants.values.length).to eql(1)
        # Change this based on source model.  I hard coded my local mvir1 to have 3 rather than 80K for faster development.
        expect(changeset.matrix_constants.values.first.length).to eql(3)
        changeset.matrix_constants.keys.each do |k|
          changeset.matrix_constants[k] = ['this', 'is', 'matrix', 'constant']
        end

        sync_workflow = workflow.plan_synchronization(changeset, test_project, 'timepoint')
        # sync_workflow.update_block = Proc.new { |a| pp a }

        sync_workflow.execute_planned!

        resulting_models = magma_client.retrieve(Etna::Clients::Magma::RetrievalRequest.new(project_name: test_project)).models

        expect(resulting_models.model_keys.sort).to eql([
            "admission_lab", "analyte", "clinical_lab", "comorbidity",
            "cytof", "cytof_pool", "cytokine", "document", "immunoassay",
            "olink", "patient", "project",
            "rna_seq", "rna_seq_plate", "sc_rna_seq", "sc_rna_seq_pool", "symptom",
            "timepoint", "treatment"
        ])

        # Ensure, however, that linked and non parent models do not have attributes filled out.  It only includes skeleton, identifier, and parent.
        expect(resulting_models.model('symptom').template.attributes.attribute_keys.sort).to eql(['created_at', 'id', 'patient', 'updated_at'])
        expect(resulting_models.model('cytof_pool').template.attributes.attribute_keys.sort).to eql(["created_at", "cytof", "pool_name", "project", "updated_at"])

        # Ensure, however, that the core model and its tree fill out
        expect(resulting_models.model('timepoint').template.attributes.attribute_keys.length).to eql(40)
        expect(resulting_models.model('rna_seq').template.attributes.attribute_keys.length).to eql(72)
        expect(resulting_models.model('rna_seq').template.attributes.attribute('gene_tpm').options).to eql(["this", "is", "matrix", "constant"])
      end
    end

    it 'can sync a full project from point A to point B' do
      configure_etna_yml

      VCR.use_cassette('add_project_models_workflow-full-project.e2e') do
        # Change this name when re-recording the cassette file to ensure a new project is synced
        test_project = "test_add_project_models_workflow_full_tzh"
        magma_client = Etna::Clients::Magma.new(
            host: 'https://magma.development.local',
            token: ENV['TOKEN'] || TEST_TOKEN, ignore_ssl: true,
        )

        magma_client.update_model(Etna::Clients::Magma::UpdateModelRequest.new(
            project_name: test_project,
            actions: [Etna::Clients::Magma::AddProjectAction.new]))

        # Ensure this is a brand new project we are testing with.
        expect(magma_client.retrieve(Etna::Clients::Magma::RetrievalRequest.new(project_name: test_project)).models.model_keys).to eql(['project'])

        workflow = Etna::Clients::Magma::AddProjectModelsWorkflow.new(magma_client: magma_client)

        io = StringIO.new
        workflow.write_models_template_csv('mvir1', io: io)
        io.rewind

        changeset = workflow.prepare_changeset_from_csv(io: io) do |err|
          raise err
        end

        # Change this based on source model
        expect(changeset.matrix_constants.values.length).to eql(1)
        # Change this based on source model.  I hard coded my local mvir1 to have 3 rather than 80K for faster development.
        expect(changeset.matrix_constants.values.first.length).to eql(3)
        changeset.matrix_constants.keys.each do |k|
          changeset.matrix_constants[k] = ['this', 'is', 'matrix', 'constant']
        end

        sync_workflow = workflow.plan_synchronization(changeset, test_project)
        # sync_workflow.update_block = Proc.new { |a| pp a }

        sync_workflow.execute_planned!

        resulting_models = magma_client.retrieve(Etna::Clients::Magma::RetrievalRequest.new(project_name: test_project)).models
        expect(resulting_models.model_keys.sort).to eql([
          "admission_lab", "analyte", "clinical_lab", "comorbidity",
          "cytof", "cytof_pool", "cytokine", "document", "immunoassay",
          "olink", "patient", "project",
          "rna_seq", "rna_seq_plate", "sc_rna_seq", "sc_rna_seq_pool", "symptom",
          "timepoint", "treatment"
        ])

        # Ensure, unlike the above partial sync test, that the entire tree and all attributes are synced.
        expect(resulting_models.model('symptom').template.attributes.attribute_keys.sort).to eql(["bleeding_site", "created_at", "id", "name", "other_name", "patient", "present", "symptom_date", "symptom_reported_by", "updated_at"])
        expect(resulting_models.model('cytof_pool').template.attributes.attribute_keys.length).to eql(17)
        expect(resulting_models.model('timepoint').template.attributes.attribute_keys.length).to eql(40)
        expect(resulting_models.model('rna_seq').template.attributes.attribute_keys.length).to eql(72)
        expect(resulting_models.model('rna_seq').template.attributes.attribute('gene_tpm').options).to eql(["this", "is", "matrix", "constant"])
      end
    end
  end
end
