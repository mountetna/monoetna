describe Vulcan::Snakemake::Utils do
  include Rack::Test::Methods
  let(:dag) { ['job1', 'job2', 'job3', 'job4', 'job5'] }

  it 'returns jobs after the last matched job' do
    jobs = ['job2', 'job3']
    result = Vulcan::Snakemake::Utils.find_affected_downstream_jobs(dag, jobs)
    expect(result).to eq(['job4', 'job5'])
  end

  it 'returns an empty array if the last job is the last element' do
    jobs = ['job5']
    result = Vulcan::Snakemake::Utils.find_affected_downstream_jobs(dag, jobs)
    expect(result).to eq([])
  end

  it 'raises an error indicating no matching elements found' do
    jobs = ['job6']
    expect { Vulcan::Snakemake::Utils.find_affected_downstream_jobs(dag, jobs) }.to raise_error("Cannot find any matching jobs in the dag")
  end
end