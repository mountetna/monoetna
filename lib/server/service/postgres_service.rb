module PostgresService
  def self.create_new_file!(params)
    FileModel::File.create(
      original_name: params['original_name'],
      file_name: params['file_name'],
      file_size: params['file_size'],
      group_name: params['group_name'],
      project_name: params['project_name'],
      project_name_full: params['project_name_full'],
      start_upload: Time.at(params['start_timestamp'].to_i),
      finish_upload: params['finish_upload'],
      upload_by: params['user_email'],
      hashing_algorithm: params['hashing_algorithm'],
      hash: params['hash']
    )
  end

  def self.create_new_upload!(params, file_id)
    FileModel::Upload.create(
      current_byte_position: params['current_byte_position'],
      current_blob_size: params['current_blob_size'],
      next_blob_size: params['next_blob_size'],
      next_blob_hash: params['next_blob_hash'],
      status: params['status'],
      file_id: file_id
    )
  end
end
