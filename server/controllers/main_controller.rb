class MainController < BasicController

  def run()

    set_user()
    raise_err(:BAD_REQ, 2, __method__) if !@user.valid?()

    return send(@action).to_json()
  end

  def retrieve_files()

    project_names = extract_project_names()
    file_list = pull_file_metadata(project_names)
    file_list = pull_upload_metadata(file_list)
    { :success=> true, :file_list=> file_list }
  end

  private
  def extract_project_names()

    @user.permissions.map do |permission|

      permission['project_name']
    end
  end

  def pull_file_metadata(project_names)

    file_metadata = []
    project_names.each do |project_name|

      project_files = PostgresService::get_files_by_project_name(project_name)
      if project_files then file_metadata.concat(project_files) end
    end
    return file_metadata
  end

  def pull_upload_metadata(file_list)

    file_list.map do |file|

      if file[:upload_by] == @user.email()

        upload = FileModel::Upload[:file_id=> file[:file_id]]
        if upload then file.concat(upload) end
      end
      file
    end
  end
end
