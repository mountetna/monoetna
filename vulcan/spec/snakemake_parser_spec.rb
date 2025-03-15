describe Vulcan::Snakemake::TargetParser do
  include Rack::Test::Methods

  let(:core_snakefile) { File.read('/app/spec/fixtures/v2/snakemake-repo/rules/core.smk') }
  let(:summary_snakefile) { File.read('/app/spec/fixtures/v2/snakemake-repo/rules/summary-ui.smk') }
  let(:config_yaml) { YAML.load_file('/app/spec/fixtures/v2/snakemake-repo/default-config.json') }
  let(:parsed_core) do
    parser = Vulcan::Snakemake::TargetParser.new(core_snakefile, config_yaml)
    parser.parse
  end
  let(:parsed_summary) do
    parser = Vulcan::Snakemake::TargetParser.new(summary_snakefile, config_yaml)
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
        expect(parsed_core["output/arithmetic.txt"]["inputs"]).to contain_exactly("output/count_poem.txt", "output/count_poem_2.txt")
        expect(parsed_core["output/check.txt"]["inputs"]).to contain_exactly("output/arithmetic.txt")
  end

    it 'correctly maps the targets with the params' do
        expect(parsed_core["output/count_poem.txt"]["params"]).to contain_exactly("count_bytes", "count_chars")
        expect(parsed_core["output/count_poem_2.txt"]["params"]).to contain_exactly("count_bytes", "count_chars")
        expect(parsed_core["output/arithmetic.txt"]["params"]).to contain_exactly("add", "add_and_multiply_by")
    end

  end

  context 'summary snakefile' do

    it 'does not include targets from ui rules' do
      require 'pry'; binding.pry
      expect(parsed_summary.keys).not_to include("output/ui_job_one.txt")
      expect(parsed_summary.keys).not_to include("output/ui_job_two.txt")
    end

    it 'correctly maps the targets with the inputs' do
      expect(parsed_summary["output/summary.txt"]["inputs"]).to contain_exactly(
        "output/count_poem.txt", "output/count_poem_2.txt",
        "output/arithmetic.txt", "output/check.txt",
        "output/ui_job_one.txt", "output/ui_job_two.txt"
      )
    end
  end
end
