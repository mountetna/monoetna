class DownloadController < Metis::Controller
  # This is the endpoint that allows you to make a download.
  # You may call this with a token
  def authorize
  end

  def download
    bucket = require_bucket

    file = Metis::File.from_path(bucket, @params[:file_path])

    raise Etna::Error.new('File not found', 404) unless file && file.has_data?

    return [
      200,
      { 'X-Sendfile' => file.data_block.location },
      [ '' ]
    ]
  end
end
