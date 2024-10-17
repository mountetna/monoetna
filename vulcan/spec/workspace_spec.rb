describe Vulcan::Workspace do
    let(:core_snakefile) { File.read('/app/spec/fixtures/v2/snakemake-repo/rules/core.smk') }
    let(:summary_snakefile) { File.read('/app/spec/fixtures/v2/snakemake-repo/rules/summary-ui.smk') }
    let(:config_yaml) { YAML.load_file('/app/spec/fixtures/v2/snakemake-repo/config.yaml') }
  
    let(:target_mapping) do
        core_parser = Vulcan::Snakemake::TargetParser.new(core_snakefile, config_yaml)
        summary_parser = Vulcan::Snakemake::TargetParser.new(summary_snakefile, config_yaml)
        core_parser.parse.merge(summary_parser.parse)
    end

    let(:workflow) do
        Vulcan::WorkflowV2.create(
            project_name: 'Test Project',
            name: 'Test Workflow',
            repo_remote_url: "test@remote-repo.com",
            created_at: Time.now,
            updated_at: Time.now
        )
    end

    let(:workspace) do
        Vulcan::Workspace.create(
            target_mapping: target_mapping,
            workflow_id: workflow.id, 
            name: 'Test Workspace',
            user_email: 'test@example.com',
            path: '/app/spec/workspaces/test_workspace',
            tags: [],
            dag: [],
            created_at: Time.now,
            updated_at: Time.now
        )
    end

    it 'correctly finds the first two targets' do
        provided_params =  ["count_bytes", "count_chars"]
        available_files = ["output/poem.txt", "output/poem_2.txt"]
        targets = workspace.find_buildable_targets(provided_params, available_files)
        expect(targets).to eq((["output/count_poem.txt", "output/count_poem_2.txt"]))

    end
    
    it 'correctly finds targets that depend on the previous targets' do
        # Here "output/count_poem.txt", "output/count_poem_2.txt" are buildable, but 
        # once they are built they contribute to the list of available files 
        # and "output/arithmetic.txt", "output/check.txt" become buildable.
        provided_params =  ["count_bytes", "count_chars", "add", "add_and_multiply_by" ]
        available_files = ["output/poem.txt", "output/poem_2.txt"]
        targets = workspace.find_buildable_targets(provided_params, available_files)
        expect(targets).to eq((["output/count_poem.txt", "output/count_poem_2.txt", "output/arithmetic.txt", "output/check.txt"]))
    end

    it 'correctly finds all targets with all params' do
        provided_params =  ["count_bytes", "count_chars", "add", "add_and_multiply_by" ]
        available_files = ["output/poem.txt", "output/poem_2.txt", "output/ui_job_one.txt", "output/ui_job_two.txt"]
        targets = workspace.find_buildable_targets(provided_params, available_files)
        expect(targets).to eq((["output/count_poem.txt", "output/count_poem_2.txt", "output/arithmetic.txt", "output/check.txt", "output/summary.txt"]))
    end

end
