describe Vulcan::Snakemake::TargetParser do
  include Rack::Test::Methods

  let(:core_snakefile) { File.read('/app/spec/fixtures/v2/snakemake-repo/rules/core.smk') }
  let(:summary_snakefile) { File.read('/app/spec/fixtures/v2/snakemake-repo/rules/summary-ui.smk') }
  let(:config_yaml) { JSON.parse(File.read('/app/spec/fixtures/v2/snakemake-repo/default-config.json')) }
  let(:parsed_core) do
    parser = Vulcan::Snakemake::TargetParser.new(core_snakefile, config_yaml)
    parser.parse
  end
  let(:parsed_summary) do
    parser = Vulcan::Snakemake::TargetParser.new(summary_snakefile, config_yaml)
    parser.parse
  end

  let(:sc_viz_snakefile) { File.read('/app/spec/fixtures/v2/other_snakefiles/sc_viz.smk') }
  let(:sc_viz_config) { JSON.parse(File.read('/app/spec/fixtures/v2/other_snakefiles/default_config.json')) }
  let(:parsed_sc_viz) do
    parser = Vulcan::Snakemake::TargetParser.new(sc_viz_snakefile, sc_viz_config)
    parser.parse
  end

  context 'core snakefile' do

    it 'correctly replaces the config variables in the snakefile' do
      expect(parsed_core["output/count_poem.txt"]["inputs"]).to contain_exactly("output/poem.txt", "output/poem_2.txt")
      expect(parsed_core["output/count_poem_2.txt"]["inputs"]).to contain_exactly("output/poem.txt", "output/poem_2.txt")
    end

    it 'correctly parses the targets' do
        expect(parsed_core.keys).to contain_exactly(
        "output/count_poem.txt",
        "output/count_poem_2.txt",
        "output/arithmetic.txt",
        "output/check.txt"
        )
    end

    it 'correctly maps the targets with the inputs' do
        expect(parsed_core["output/count_poem.txt"]["inputs"]).to contain_exactly("output/poem.txt", "output/poem_2.txt")
        expect(parsed_core["output/count_poem_2.txt"]["inputs"]).to contain_exactly("output/poem.txt", "output/poem_2.txt")
        expect(parsed_core["output/arithmetic.txt"]["inputs"]).to contain_exactly("output/count_poem.txt", "output/count_poem_2.txt", "resources/number_to_add.txt")
        expect(parsed_core["output/check.txt"]["inputs"]).to contain_exactly("output/arithmetic.txt")
  end

    it 'correctly maps the targets with the params' do
        expect(parsed_core["output/count_poem.txt"]["params"]).to contain_exactly("count_bytes", "count_chars")
        expect(parsed_core["output/count_poem_2.txt"]["params"]).to contain_exactly("count_bytes", "count_chars")
        expect(parsed_core["output/arithmetic.txt"]["params"]).to contain_exactly("add", "add_and_multiply_by")
    end

    it 'correctly maps the targets with the rule names' do
        expect(parsed_core["output/count_poem.txt"]["rule"]).to eq("count")
        expect(parsed_core["output/count_poem_2.txt"]["rule"]).to eq("count")
        expect(parsed_core["output/arithmetic.txt"]["rule"]).to eq("arithmetic")
        expect(parsed_core["output/check.txt"]["rule"]).to eq("checker")
    end

  end

  context 'sc_viz snakefile' do

    it 'correctly ignores any directive that is not input, output, or params' do
      expect(parsed_sc_viz.keys).not_to include(
        "/dscolab/vulcan/containers/archimedes-r.sif",
        "scripts/get_dataset_and_summarize.R", 
        "scripts/make_dittoSeq_plot.R"
      )
    end

    it 'correctly makess sure each target only has input and params' do
      parsed_sc_viz.each do |target, properties|
        expect(properties.keys).to contain_exactly("inputs", "params", "rule")
      end
    end

    it 'uses the config variable names as params instead of named params' do
      expect(parsed_sc_viz["output/scdata.Rds"]["params"]).to contain_exactly("dataset_name")
    end

    it 'correctly maps the targets with the rule names' do
      expect(parsed_sc_viz["output/scdata.Rds"]["rule"]).to eq("get_dataset_and_summarize")
      expect(parsed_sc_viz["output/plotting_options.json"]["rule"]).to eq("get_dataset_and_summarize")
      expect(parsed_sc_viz["output/discrete_metadata_summary.json"]["rule"]).to eq("get_dataset_and_summarize")
      expect(parsed_sc_viz["output/plot_setup.json"]["rule"]).to eq("plot_setup_ui")
      expect(parsed_sc_viz["output/plot.out"]["rule"]).to eq("make_plot")
      expect(parsed_sc_viz["output/thumbnail.png"]["rule"]).to eq("make_plot")
      expect(parsed_sc_viz["output/plot.Rds"]["rule"]).to eq("make_plot")
    end

  end

  context 'summary snakefile' do

    it 'includes params ui for ui job files' do
      expect(parsed_summary["output/ui_job_one.txt"]["params"]).to include("ui")
      expect(parsed_summary["output/ui_job_two.txt"]["params"]).to include("ui")
      expect(parsed_summary["output/ui_summary.txt"]["params"]).to include("ui")
    end

    it 'correctly maps the targets with the inputs' do
      expect(parsed_summary["output/summary.txt"]["inputs"]).to contain_exactly(
        "output/count_poem.txt", "output/count_poem_2.txt",
        "output/arithmetic.txt", "output/check.txt",
        "output/ui_job_one.txt", "output/ui_job_two.txt"
      )
      expect(parsed_summary["output/ui_job_one.txt"]["inputs"]).to contain_exactly("output/check.txt")
      expect(parsed_summary["output/ui_job_two.txt"]["inputs"]).to be_empty
      expect(parsed_summary["output/ui_summary.txt"]["inputs"]).to contain_exactly("output/ui_job_one.txt", "output/ui_job_two.txt","output/summary.txt")
      expect(parsed_summary["output/final.txt"]["inputs"]).to contain_exactly("output/ui_summary.txt")
    end

    it 'correctly maps the targets with the rule names' do
      expect(parsed_summary["output/ui_job_one.txt"]["rule"]).to eq("ui_job_one")
      expect(parsed_summary["output/ui_job_two.txt"]["rule"]).to eq("ui_job_two")
      expect(parsed_summary["output/ui_summary.txt"]["rule"]).to eq("ui_summary")
      expect(parsed_summary["output/summary.txt"]["rule"]).to eq("summary")
      expect(parsed_summary["output/final.txt"]["rule"]).to eq("final")
    end
  end
end
