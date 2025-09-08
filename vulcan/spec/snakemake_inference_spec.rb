describe Vulcan::Snakemake::Inference do
  let(:core_snakefile) { File.read('/app/spec/fixtures/v2/snakemake-repo/rules/core.smk') }
  let(:summary_snakefile) { File.read('/app/spec/fixtures/v2/snakemake-repo/rules/summary-ui.smk') }
  let(:config_yaml) { YAML.load_file('/app/spec/fixtures/v2/snakemake-repo/default-config.json') }

  let(:target_mapping) do
    core_parser = Vulcan::Snakemake::TargetParser.new(core_snakefile, config_yaml)
    summary_parser = Vulcan::Snakemake::TargetParser.new(summary_snakefile, config_yaml)
    core_parser.parse.merge(summary_parser.parse)
  end


  context 'find_buildable_targets' do
    it 'correctly finds the first two targets' do
      provided_params = ["count_bytes", "count_chars"]
      available_files = ["output/poem.txt", "output/poem_2.txt"]
      targets = Vulcan::Snakemake::Inference.find_buildable_targets(target_mapping, provided_params, available_files)
      expect(targets).to eq(["output/count_poem.txt", "output/count_poem_2.txt"])
    end

    it 'correctly finds targets that depend on the previous targets' do
      # Here "output/count_poem.txt", "output/count_poem_2.txt" are buildable, but 
      # once they are built they contribute to the list of available files 
      # and "output/arithmetic.txt", "output/check.txt" become buildable.
      provided_params = ["count_bytes", "count_chars", "add", "add_and_multiply_by"]
      available_files = ["output/poem.txt", "output/poem_2.txt", "resources/number_to_add.txt"]
      targets = Vulcan::Snakemake::Inference.find_buildable_targets(target_mapping, provided_params, available_files)
      expect(targets).to eq(["output/count_poem.txt", "output/count_poem_2.txt", "output/arithmetic.txt", "output/check.txt"])
    end

    it 'correctly finds all targets with all params' do
      provided_params = ["count_bytes", "count_chars", "add", "add_and_multiply_by", "ui"]
      available_files = ["output/poem.txt", "output/poem_2.txt", "resources/number_to_add.txt", "output/ui_job_one.txt", "output/ui_job_two.txt"]
      targets = Vulcan::Snakemake::Inference.find_buildable_targets(target_mapping, provided_params, available_files)
      expect(targets).to eq(["output/count_poem.txt", "output/count_poem_2.txt", "output/arithmetic.txt", "output/check.txt", "output/summary.txt"])
    end
  end 

  context 'find_ui_targets' do
    it 'correctly finds the ui targets' do
      result = Vulcan::Snakemake::Inference.ui_targets(target_mapping)
      expect(result).to eq(["output/ui_job_one.txt", "output/ui_job_two.txt", "output/ui_summary.txt"])
    end
  end

  context 'file graph' do
    it 'correctly builds the adjacency list' do
      result = Vulcan::Snakemake::Inference.file_graph(target_mapping)
      expect(result).to match(
          a_hash_including(
              "output/poem.txt"=>["output/count_poem.txt", "output/count_poem_2.txt"],
              "output/poem_2.txt"=>["output/count_poem.txt", "output/count_poem_2.txt"],
              "output/count_poem.txt"=>["output/arithmetic.txt", "output/summary.txt"],
              "output/count_poem_2.txt"=>["output/arithmetic.txt", "output/summary.txt"],
              "output/arithmetic.txt"=>["output/check.txt", "output/summary.txt"],
              "output/check.txt"=>["output/ui_job_one.txt", "output/summary.txt"],
              "output/ui_job_one.txt"=>["output/summary.txt", "output/ui_summary.txt"],
              "output/ui_job_two.txt"=>["output/summary.txt", "output/ui_summary.txt"],
              "output/summary.txt"=>["output/ui_summary.txt"],
              "output/ui_summary.txt"=>["output/final.txt"]
          )
      )
    end
  end

  context 'downstream_nodes' do
    let(:file_graph) { Vulcan::Snakemake::Inference.file_graph(target_mapping) }
  
    it 'correctly finds downstream nodes' do
      result = Vulcan::Snakemake::Inference.downstream_nodes(file_graph, ["output/summary.txt"])
      expect(result).to eq(Set.new(["output/ui_summary.txt", "output/final.txt"]))
    end

    it 'correctly finds downstream nodes when there are multiple targets' do
      result = Vulcan::Snakemake::Inference.downstream_nodes(file_graph, ["output/arithmetic.txt", "output/ui_job_two.txt"])
      expect(result).to eq(Set.new(["output/check.txt", "output/summary.txt", "output/ui_summary.txt", "output/ui_job_one.txt", "output/final.txt"]))
    end

    it 'correctly ignores downstream nodes on the same level' do
      result = Vulcan::Snakemake::Inference.downstream_nodes(file_graph, ["output/ui_job_one.txt"])
      # ui_job_two.txt is on the same level as ui_job_one.txt
      expect(result).to eq(Set.new(["output/summary.txt", "output/ui_summary.txt", "output/final.txt"]))
    end
  end

  end
end