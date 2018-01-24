class DownloadController < Metis::Controller
  # This is the endpoint that allows you to make a download.
  # You may call this with a token
  def authorize
  end

  def download
    file = Metis::File.where(
      project_name: @params[:project_name],
      file_name: @params[:file_name]
    ).first

    return failure(404, 'File not found') unless file && file.has_data?

    return [
      200,
      { 'X-Sendfile' => file.location },
      ''
    ]
  end
end
