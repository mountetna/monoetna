describe Vulcan::Snakemake::TargetParser do
  include Rack::Test::Methods

  let(:core_snakefile) { File.read('/app/spec/fixtures/v2/snakemake-repo/rules/core.smk') }
  let(:summary_snakefile) { File.read('/app/spec/fixtures/v2/snakemake-repo/rules/summary-ui.smk') }
  let(:config_yaml) { YAML.load_file('/app/spec/fixtures/v2/snakemake-repo/config.yaml') }
  let(:parsed) do
    parser = Vulcan::Snakemake::TargetParser.new(core_snakefile, config_yaml)
    parser.parse
  end

  it 'correctly replaces the config variables in the snakefile' do
    expect(parsed["output/count_poem.txt"][:inputs]).to contain_exactly("output/poem.txt", "output/poem_2.txt")
    expect(parsed["output/count_poem_2.txt"][:inputs]).to contain_exactly("output/poem.txt", "output/poem_2.txt")
  end

  it 'correctly parses the targets' do
    expected_keys = [
      "output/count_poem.txt",
      "output/count_poem_2.txt",
      "output/arithmetic.txt",
      "output/check.txt"
    ]
    expect(parsed.keys).to match_array(expected_keys)
  end


  it 'correctly maps the targets with the inputs' do
    expect(parsed["output/count_poem.txt"][:inputs]).to contain_exactly("output/poem.txt", "output/poem_2.txt")
    expect(parsed["output/count_poem_2.txt"][:inputs]).to contain_exactly("output/poem.txt", "output/poem_2.txt")
    expect(parsed["output/arithmetic.txt"][:inputs]).to contain_exactly("output/count_poem.txt", "output/count_poem_2.txt")
    expect(parsed["output/check.txt"][:inputs]).to contain_exactly("output/arithmetic.txt")
  end

  it 'correctly maps the targets with the params' do
    expect(parsed["output/count_poem.txt"][:params]).to contain_exactly("count_bytes", "count_chars")
    expect(parsed["output/count_poem_2.txt"][:params]).to contain_exactly("count_bytes", "count_chars")
    expect(parsed["output/arithmetic.txt"][:params]).to contain_exactly("add", "add_and_multiply_by")
    expect(parsed["output/check.txt"][:params]).to be_empty
  end

  it 'correctly does not include targets from ui rules' do
    parser = Vulcan::Snakemake::TargetParser.new(summary_snakefile, config_yaml)
    parsed_summary = parser.parse
    expect(parsed_summary.keys).not_to include("output/checker-ui.txt")
  end

end
