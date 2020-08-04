require 'mimemagic'
require 'vips'

class DownloadController < Metis::Controller
  # This is the endpoint that allows you to make a download.
  # You may call this with a token
  def authorize
  end

  def download
    bucket = require_bucket

    file = Metis::File.from_path(bucket, @params[:file_path])

    raise Etna::Error.new('File not found', 404) unless file && file.has_data?

    puts @params
    if @params.fetch(:thumbnail, nil)
      puts file.data_block.location
      thumbnail = Vips::Image.thumbnail(file.data_block.location, 240)
      puts thumbnail
      puts thumbnail.filename
      puts MimeMagic.by_path(file.file_name)
      return [
        200,
        { 'Content-Type' => MimeMagic.by_path(file.file_name) },
        [ thumbnail.tiffsave_buffer ]
      ]
    end

    return [
      200,
      { 'X-Sendfile' => file.data_block.location },
      [ '' ]
    ]
  end
end
