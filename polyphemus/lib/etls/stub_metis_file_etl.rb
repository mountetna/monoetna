require_relative '../metis_file_etl'

class Polyphemus::StubMetisFileEtl < Polyphemus::MetisFileEtl
  def initialize
    super(project_bucket_pairs: [])
  end

  def process(cursor, files)
  end
end