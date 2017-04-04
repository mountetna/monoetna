class MainController < BasicController

  def run()

    set_user()
    raise_err(:BAD_REQ, 2, __method__) if !@user.valid?()

    return send(@action).to_json()
  end

  def retrieve_files()

    project_names = extract_project_names()
    file_list = pull_file_metadata(project_names)
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

      #puts project_name
    end

    return []
  end
end