class FileModel

  def initialize(basic_metadata = nil, redis_service = nil)

    @metadata = basic_metadata
    @redis_service = redis_service

    @file_exists = false
    @metadata_exists = false

    @db_index = nil
    @directory = nil
    @file_name = nil
    @full_path = nil
    @metadata_key = nil

    set_basic_file_data()
    check_existance()
    #verify_data()
  end

  def file_exists?()

    @file_exists
  end

  def metadata_exists?()

    @metadata_exists
  end

  def db_index()

    @db_index
  end

  def directory()

    @directory
  end

  def file_name()

    @file_name
  end

  def group_id()

    @group_id
  end

  def project_id()

    @project_id
  end

  # Generate the key used to access the file's metadata in Redis.
  def metadata_key()

    @metadata_key
  end

  private

  # Set the file model based on basic metadata.
  def set_basic_file_data()

    # Bail out of setting up this model if the input params are not correct.
    if !has_required_params?() then return end

    @group_id = @metadata['group_id']
    @project_id = @metadata['project_id']
    @directory = Conf::ROOT_DIR+'/'+@group_id.to_s+'/'+@project_id.to_s()
    @file_name = @metadata['file_name']
    @full_path = @directory+'/'+@file_name
  end

  # Check that the incoming basic metadata has the correct parameters.
  def has_required_params?()

    has_params = true
    ['project_name', 'project_id', 'file_name', 'group_id'].each do |param|

      if !@metadata.key?(param) then has_params = false end
    end
    return has_params
  end

  # Check for the exsistance of a file on disk and metadata in Redis.
  def check_existance()

    # Check for the file on disk.
    if File.file?(@full_path) then @file_exists = true end

    # Check for the file metadata.
    partial_key = '*.'+@file_name+'.'+@group_id.to_s+'.'+@project_id.to_s
    redis_keys = @redis_service.retrieve_file_key(partial_key)
    if redis_keys.length > 0

      @metadata_key = redis_keys[0]
      @metadata_exists = true
    end
  end

  # Check that the incoming basic metadata matches with the Redis metadata.
  #def verify_data()

  #end
end