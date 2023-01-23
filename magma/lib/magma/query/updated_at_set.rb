require_relative "metis_metadata"

class UpdatedAtSet < MetisMetadata
  def [](file_path)
    super

    @requested_files[file_path] ?
      Magma::Answer.new(@requested_files[file_path][:updated_at]) :
      Magma::NilAnswer.new
  end
end
