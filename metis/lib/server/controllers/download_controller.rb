class DownloadController < Metis::Controller
  # This is the endpoint that allows you to make a download.
  # You may call this with a token
  def authorize
    bucket = require_bucket
    success_json(download_url: Metis::File.download_url(@request, bucket.project_name, bucket.name, @params[:file_path]))
  end

  def download
    bucket = require_bucket

    file = Metis::File.from_path(bucket, @params[:file_path])

    raise Etna::Error.new('File not found', 404) unless file && file.has_data?

    if @params.fetch(:thumbnail, nil) && file.data_block.thumbnail_in_cache?
      return [
        200,
        {},
        [ file.data_block.thumbnail ]
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
  rescue Metis::ThumbnailError => e
    Metis.instance.logger.log_error(e)
    raise Etna::Error.new("Error generating thumbnail, #{e.message}.", 422)
  end
end
