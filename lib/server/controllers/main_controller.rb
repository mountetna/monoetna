class MainController < Metis::Controller
  def retrieve_files
    files = Metis::Files.where(project_name: @user.projects).all
    uploads = files.select { |f| f.uploader == @user.email && f.upload }.map(&:upload)
    success({ files: files, uploads: uploads }.to_json, 'application/json')
  end
end
