class MainController < BasicController

  def run()

    set_user()
    if !@user.valid?() then raise_err(:BAD_REQ, 2, __method__) end

    return send(@action).to_json()
  end

  def retrieve_files()

    metadata_ids = extract_metadata_ids()
    file_list = pull_file_metadata(metadata_ids, @user.id)
    { :success=> true, :file_list=> file_list }
  end

  private
  def extract_metadata_ids()

    @user.permissions.map do |permission|

      permission['group_id'].to_s()+'.'+permission['project_id'].to_s()
    end
  end

  def pull_file_metadata(metadata_ids, user_id)

    file_metadata_keys = []
    file_metadata = []

    metadata_ids.each do |metadata_id|

      redis_keys = @redis_service.retrieve_file_key('*'+metadata_id)
      if redis_keys != nil then file_metadata_keys.push(*redis_keys) end
    end

    file_metadata_keys.each do |file_metadata_key|

      redis_metadata = @redis_service.retrieve_file_status(file_metadata_key)
      if (redis_metadata != nil && redis_metadata.key?('finish_timestamp')) ||
         (redis_metadata['user_id'].to_s == user_id.to_s)

        file_metadata.push(redis_metadata)
      end
    end

    return file_metadata
  end
end