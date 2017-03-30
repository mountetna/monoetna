class UserModel

  def initialize(user_data = nil)

    @user_data = user_data
    @valid = false

    @email = nil
    @first_name = nil
    @last_name = nil
    @user_id = nil
    @token = nil
    @permissions = nil

    validate_user_data()
    set_user_data()
  end

  def valid?()

    @valid
  end

  def email()

    @email
  end

  def first_name()

    @first_name
  end

  def last_name()

    @last_name
  end

  def id()

    @user_id
  end

  def token()

    @token
  end

  def permissions()

    @permissions
  end

  def project_editor?(project_name)

    editor = false
    @permissions.each do |perm|

      if perm['project_name'].to_s == project_name.to_s &&
         perm['role'] == 'editor'

        editor = true
      end
    end
    return editor
  end

  def project_admin?(project_name)

    admin = false
    @permissions.each do |perm|

      if perm['project_name'].to_s == project_name.to_s &&
         perm['role'] == 'administrator'

        admin = true
      end
    end
    return admin
  end

  private
  def validate_user_data()

    m = __method__

    # Check for the appropriate data fields
    if !@user_data.key?('user_info') then raise BasicError.new(:BAD_REQ,0,m) end
    user_info = @user_data['user_info']
    fields = ['email','first_name','last_name','user_id','token','permissions']
    fields.each do |field|

      if !user_info.key?(field) then raise BasicError.new(:BAD_REQ, 0, m) end
    end

    @valid = true
  end

  def set_user_data()

    unless @valid then return end

    @email = @user_data['user_info']['email']
    @first_name = @user_data['user_info']['first_name']
    @last_name = @user_data['user_info']['last_name']
    @user_id = @user_data['user_info']['user_id']
    @token = @user_data['user_info']['token']
    @permissions = @user_data['user_info']['permissions']
  end
end