require_relative '../../thumbnail_cache'

class DownloadController < Metis::Controller
  # This is the endpoint that allows you to make a download.
  # You may call this with a token
  def authorize
  end

  def download
    bucket = require_bucket

    file = Metis::File.from_path(bucket, @params[:file_path])

    raise Etna::Error.new('File not found', 404) unless file && file.has_data?

    @cache = Metis::ThumbnailCache.new

    if @params.fetch(:thumbnail, nil) && @cache.thumbnail_in_cache?(file)
      return [
        200,
        { 'Content-Type' => @cache.mimetype(file) },
        [ @cache.thumbnail(file) ]
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
