require_relative "metis_metadata"

class Md5Set < MetisMetadata
  def [](file_path)
    super

    @requested_files[file_path] ?
      Magma::Answer.new(@requested_files[file_path][:file_hash]) :
      Magma::NilAnswer.new
  end
end
