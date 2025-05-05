require 'webmock/rspec'
require 'json'
require_relative '../../../lib/etna/clients/metis'

describe 'Metis Client class' do
  let(:test_class) { Etna::Clients::Metis.new(token: TEST_TOKEN, host: METIS_HOST) }

  before(:each) do
    stub_metis_setup
  end

  it 'can fetch byte counts' do
    stub_metis_get_byte_count

    test_class.get_byte_count_by_project(Etna::Clients::Metis::GetFileCountByProjectRequest.new)
    expect(WebMock).to have_requested(:get, %r!#{METIS_HOST}/api/stats/bytes!)
  end

  it 'can fetch file counts with projects specified' do
    project_names = ['test']

    stub_metis_get_file_count(project_names)

    query = project_names.map do |name|
      "projects%5B%5D=#{name}"
    end.join('&')

    test_class.get_file_count_by_project(Etna::Clients::Metis::GetFileCountByProjectRequest.new(project_names: project_names))
    expect(WebMock).to have_requested(:get, %r!#{METIS_HOST}/api/stats/files\?#{query}!)
  end

  it 'can fetch byte counts' do
    stub_metis_get_byte_count

    test_class.get_byte_count_by_project(Etna::Clients::Metis::GetByteCountByProjectRequest.new)
    expect(WebMock).to have_requested(:get, %r!#{METIS_HOST}/api/stats/bytes!)
  end

  it 'can fetch byte counts with projects specified' do
    project_names = ['test']

    stub_metis_get_byte_count(project_names)

    query = project_names.map do |name|
      "projects%5B%5D=#{name}"
    end.join('&')

    test_class.get_byte_count_by_project(Etna::Clients::Metis::GetByteCountByProjectRequest.new(project_names: project_names))
    expect(WebMock).to have_requested(:get, %r!#{METIS_HOST}/api/stats/bytes\?#{query}!)
  end
end
