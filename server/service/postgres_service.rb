module PostgresService

  def self.connect()

    db_config = {

      :adapter=> 'postgres',
      :host=> 'localhost', 
      :database=> 'metis',

      :user=> Secrets::PSQL_USER,
      :password=> Secrets::PSQL_PASS,
      :search_path=>['private']
    }

    @postgres = Sequel.connect(db_config)
  end

  def self.create_new_file!(params)

    row = {

      :original_name=> params['original_name'],
      :file_name=> params['file_name'],
      :file_size=> params['file_size'],
      :group_name=> params['group_name'],
      :project_name=> params['project_name'],
      :start_upload=> Time.at(params['start_timestamp'].to_i),
      :finish_upload=> nil,
      :upload_by=> params['user_email'],
      :hashing_algorithm=> params['hashing_algorithm'],
      :hash=> nil
    }

    FileModel::File.create(row)
  end

  def self.create_new_upload!(params, file_id)

    row = {

      :current_byte_position=> params['current_byte_position'],
      :current_blob_size=> params['current_blob_size'],
      :next_blob_size=> params['next_blob_size'],
      :next_blob_hash=> params['next_blob_hash'],
      :status=> params['status'],

      :file_id=> file_id
    }

    FileModel::Upload.create(row)
  end

  def self.get_files_by_project_name(project_name)

    file_metadata = @postgres[:files]
      .where('project_name = ?', project_name)
      .all

    if file_metadata then file_metadata.map(){ |file| file.to_hash } end
    return file_metadata
  end
end
