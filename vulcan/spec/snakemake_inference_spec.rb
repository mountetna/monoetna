describe Vulcan::Snakemake::Inference do
    let(:core_snakefile) { File.read('/app/spec/fixtures/v2/snakemake-repo/rules/core.smk') }
    let(:summary_snakefile) { File.read('/app/spec/fixtures/v2/snakemake-repo/rules/summary-ui.smk') }
    let(:config_yaml) { YAML.load_file('/app/spec/fixtures/v2/snakemake-repo/config.yaml') }

    let(:target_mapping) do
        core_parser = Vulcan::Snakemake::TargetParser.new(core_snakefile, config_yaml)
        summary_parser = Vulcan::Snakemake::TargetParser.new(summary_snakefile, config_yaml)
        core_parser.parse.merge(summary_parser.parse)
    end

    let(:dag) { ['job1', 'job2', 'job3', 'job4', 'job5'] }

    context 'find_buildable_targets' do
        it 'correctly finds the first two targets' do
        provided_params =  ["count_bytes", "count_chars"]
        available_files = ["output/poem.txt", "output/poem_2.txt"]
        targets = Vulcan::Snakemake::Inference.find_buildable_targets(target_mapping, provided_params, available_files)
        expect(targets).to eq(["output/count_poem.txt", "output/count_poem_2.txt"])
    end
    
    it 'correctly finds targets that depend on the previous targets' do
        # Here "output/count_poem.txt", "output/count_poem_2.txt" are buildable, but 
        # once they are built they contribute to the list of available files 
        # and "output/arithmetic.txt", "output/check.txt" become buildable.
        provided_params =  ["count_bytes", "count_chars", "add", "add_and_multiply_by" ]
        available_files = ["output/poem.txt", "output/poem_2.txt"]
        targets = Vulcan::Snakemake::Inference.find_buildable_targets(target_mapping, provided_params, available_files)
        expect(targets).to eq(["output/count_poem.txt", "output/count_poem_2.txt", "output/arithmetic.txt", "output/check.txt"])
    end

    it 'correctly finds all targets with all params' do
        provided_params =  ["count_bytes", "count_chars", "add", "add_and_multiply_by" ]
        available_files = ["output/poem.txt", "output/poem_2.txt", "output/ui_job_one.txt", "output/ui_job_two.txt"]
        targets = Vulcan::Snakemake::Inference.find_buildable_targets(target_mapping, provided_params, available_files)
        expect(targets).to eq(["output/count_poem.txt", "output/count_poem_2.txt", "output/arithmetic.txt", "output/check.txt", "output/summary.txt"])
    end

    end 

    context 'find_affected_downstream_jobs' do

        it 'returns jobs after the last matched job' do
            jobs = ['job2', 'job3']
            result = Vulcan::Snakemake::Inference.find_affected_downstream_jobs(dag, jobs)
            expect(result).to eq(['job4', 'job5'])
        end
          
        it 'returns an empty array if the last job is the last element' do
            jobs = ['job5']
            result = Vulcan::Snakemake::Inference.find_affected_downstream_jobs(dag, jobs)
            expect(result).to eq([])
        end
          
        it 'raises an error indicating no matching elements found' do
            jobs = ['job6']
            expect { Vulcan::Snakemake::Inference.find_affected_downstream_jobs(dag, jobs) }.to raise_error("Cannot find any matching jobs in the dag")
        end
    end
end