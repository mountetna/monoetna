class FilesController < Metis::Controller
  def index
    files = Metis::File.where(project_name: @params[:project_name]).exclude(file_hash: nil).all

    success({ files: files }.to_json, 'application/json')
  end
end
