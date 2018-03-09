class FilesController < Metis::Controller
  def index
    files = Metis::File.where(project_name: @params[:project_name])
      .exclude(file_hash: nil).all.map do |file|
      file.to_hash(@request)
    end

    success({ files: files }.to_json, 'application/json')
  end
end
