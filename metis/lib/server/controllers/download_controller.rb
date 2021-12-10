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

    if @params.fetch(:thumbnail, nil)
      thumbnail = Vips::Image.thumbnail(file.data_block.location, 240)
      mimetype = MimeMagic.by_path(file.file_name)

      case mimetype.to_s
      when /jpe?g$/
        buffer = thumbnail.jpegsave_buffer
      when 'image/tiff'
        buffer = thumbnail.tiffsave_buffer
      when 'image/png'
        buffer = thumbnail.pngsave_buffer
      else
        raise Etna::Error.new("Thumbnails not supported for mimetype #{mimetype}", 422)
      end
      return [
        200,
        { 'Content-Type' => mimetype },
        [ buffer ]
      ]
    end

    return [
      200,
      {   'X-Sendfile' => file.data_block.location,
          'Accept-Ranges' => 'bytes',
          'Content-Disposition' => "attachment; filename=#{file.file_name}",
          'Content-MD5' => file.file_hash,
      },
      [ '' ]
    ]
  rescue Vips::Error => e
    Metis.instance.logger.log_error(e)
    raise Etna::Error.new('Unknown file format -- not supported with thumbnails', 422)
  end
end
