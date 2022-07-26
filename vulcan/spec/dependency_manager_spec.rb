describe Vulcan::DependencyManager do
  let(:inputs) { {} }
  let(:figure) {
    create_figure(
      workflow_name: "test_workflow.cwl",
      dependencies: {
        "etnaagent/archimedes": "sha256:0336c5101c0a089a2bb38d2fa01c0747f4f6bd615e0326a567a256e8aa04d4b0"
      }
    )
  }
  let(:session) { Session.new_session_for('project', 'test_workflow.cwl', 'storage_key', inputs, reference_figure_id: figure.id) }
  let(:dm) { Vulcan::DependencyManager.new }

  after(:each) do
    configure_etna_yml_ignore_dependencies(false)
  end

  describe 'archimedes_run_sha' do
    it 'respects the session ignore_dependencies flag' do
      explicit_dependency = "etnaagent/archimedes@#{figure.dependencies["etnaagent/archimedes"]}"
      image_name = dm.archimedes_run_sha(session)
      expect(image_name).to eq(explicit_dependency)

      configure_etna_yml_ignore_dependencies

      ignored_image_name = dm.archimedes_run_sha(session)
      expect(ignored_image_name).not_to eq(explicit_dependency)
    end
  end

  describe 'target_image' do
    it 'respects the session ignore_dependencies flag' do
      explicit_dependency = "etnaagent/archimedes@#{figure.dependencies["etnaagent/archimedes"]}"
      image_name = dm.target_image("python", session)
      expect(image_name).to eq(explicit_dependency)

      configure_etna_yml_ignore_dependencies

      ignored_image_name = dm.target_image("python", session)
      expect(ignored_image_name).not_to eq(explicit_dependency)
    end
  end
end
