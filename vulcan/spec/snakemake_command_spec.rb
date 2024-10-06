describe Vulcan::Snakemake::CommandBuilder do
    let(:core_snakefile) { File.read('/app/spec/fixtures/v2/snakemake-repo/rules/core.smk') }
    let(:summary_snakefile) { File.read('/app/spec/fixtures/v2/snakemake-repo/rules/summary-ui.smk') }
    let(:config_yaml) { YAML.load_file('/app/spec/fixtures/v2/snakemake-repo/config.yaml') }
  
    let(:target_mapping) do
        core_parser = Vulcan::Snakemake::TargetParser.new(core_snakefile, config_yaml)
        summary_parser = Vulcan::Snakemake::TargetParser.new(summary_snakefile, config_yaml)
        core_parser.parse.merge(summary_parser.parse)
    end

    it 'correctly finds the first two targets' do
        command_builder = Vulcan::Snakemake::CommandBuilder.new
        command_builder.target_meta = {
            mapping: target_mapping,
            provided_params: ["count_bytes", "count_chars"],
            available_files: ["output/poem.txt", "output/poem_2.txt"],
        }
        targets = command_builder.determine_buildable_targets
        expect(targets).to eq((["output/count_poem.txt", "output/count_poem_2.txt"]))

    end
    
    it 'correctly finds targets that depends on the previous previous targets' do
        # Here "output/count_poem.txt", "output/count_poem_2.txt" are buildable, but 
        # once they are built they contribute to the list of available files 
        # and "output/arithmetic.txt", "output/check.txt" are become buildable.
        command_builder = Vulcan::Snakemake::CommandBuilder.new
        command_builder.target_meta = {
            mapping: target_mapping,
            provided_params: ["count_bytes", "count_chars", "add", "add_and_multiply_by" ],
            available_files: ["output/poem.txt", "output/poem_2.txt"],
        }
        targets = command_builder.determine_buildable_targets
        expect(targets).to eq((["output/count_poem.txt", "output/count_poem_2.txt", "output/arithmetic.txt", "output/check.txt"]))
    end

    it 'correctly finds all targets with all params' do
        command_builder = Vulcan::Snakemake::CommandBuilder.new
        command_builder.target_meta = {
            mapping: target_mapping,
            provided_params: ["count_bytes", "count_chars", "add", "add_and_multiply_by" ],
            available_files: ["output/poem.txt", "output/poem_2.txt", "output/ui_check.txt"],
        }
        targets = command_builder.determine_buildable_targets
        expect(targets).to eq((["output/count_poem.txt", "output/count_poem_2.txt", "output/arithmetic.txt", "output/check.txt", "output/summary.txt"]))
    end



end
