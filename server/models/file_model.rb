class FileModel

  # Define the getters and setters for this model.
  PROPERTIES = [

    # Basic info about the model that the client should always know.
    :original_name,
    :file_name,
    :file_size,                 # The expected total size of the file.
    :group_id,
    :group_name,
    :project_id,
    :project_name,

    # Basic info about the model that the server should generate.
    :directory,
    :full_path,
    :partial_file_name,         # The name of the partial file.

    # 'Complex' info about the model that the server should generate.
    :db_index,
    :metadata_key,
    :start_timestamp,
    :finish_timestamp,
    :user_id,                   # The id of the user who uploaded the file.
    :user_email,                # The email of the user who uploaded the file.
    :hashing_algorithm,

    # Info that only pertains to an file being uploaded.
    :current_byte_position,
    :current_blob_size,
    :next_blob_size,
    :next_blob_hash,
    :byte_count,                # The amount of bytes that have been uploaded.
    :status
  ]
  PROPERTIES.each { |prop| attr_accessor prop }

  def initialize(basic_metadata = {}, redis_service = nil)

    @metadata = basic_metadata
    @redis_service = redis_service

    # Info about the state of the model.
    @directory_exists = false
    @file_exists = false
    @db_metadata_exists = false
    @upload_complete = false

    # Procedures that setup and verify the model.
    set_basic_metadata()
    generate_extra_data()
    check_existance()
  end

  def directory_exists?()

    @directory_exists
  end

  def file_exists?()

    @file_exists
  end

  def db_metadata_exists?()

    @metadata_exists
  end

  def upload_complete?()

    @upload_complete
  end

  def set_file_on_system!()

    # Check that the directory exists
    if !directory_exists?() then return end

    # Check that there is no preexisting file or metadata.
    if db_metadata_exists?() || file_exists?() then return end

    initialize_metadata()
    create_partial_file()
    check_existance()
  end

  private
  def set_basic_metadata()

    @metadata.each do |key, value|

      key = key.to_sym
      if PROPERTIES.member?(key) then self.send("#{key}=", value) end
    end
  end

  def generate_extra_data()

    if self.group_name == nil ||
      self.project_name == nil ||
      self.file_name == nil then return end

    grp_nm = normalize_name(self.group_name)
    prjk_nm = normalize_name(self.project_name)

    self.directory = Conf::ROOT_DIR+'/'+grp_nm+'/'+prjk_nm
    self.full_path = self.directory+'/'+self.file_name
    self.partial_file_name = self.full_path+'.part'
  end

  # Check for the exsistance of a file on disk and metadata in Redis.
  def check_existance()

    # Check for the file's directory.
    if File.directory?(directory) then @directory_exists = true end

    # Check for the file on disk.
    if File.file?(self.full_path) then @file_exists = true end

    # Check for the file metadata.
    grp_nm = normalize_name(self.group_name)
    prjk_nm = normalize_name(self.project_name)
    redis_key = '*-'+grp_nm+'-'+prjk_nm+'-'+self.file_name
    redis_keys = @redis_service.retrieve_file_key(redis_key)
    if redis_keys.length > 0

      self.metadata_key = redis_keys[0]
      @db_metadata_exists = true
    end
  end

  # If we are here then we have had to test the HMAC first. If the HMAC 
  # is cool then all of the sensitive data here should be cool.
  def initialize_metadata()

    if @upload_complete then return end

    generate_metadata_key()

    # Modify the metadata to save.
    @metadata.delete('token')
    @metadata.delete('db_index')
    @metadata['status'] = 'active'
    @metadata['directory'] = self.directory
    @metadata['full_path'] = self.full_path
    @metadata['partial_file_name'] = self.partial_file_name

    @redis_service.set_file_status(self.metadata_key, @metadata.to_json)
  end

  # Generate the key used to access the file's metadata in Redis.
  def generate_metadata_key()
  
    items = ['db_index', 'group_name', 'project_name', 'file_name']
    items.map! do |item|
  
      if item == 'group_name' || item == 'project_name'
  
        normalize_name(@metadata[item])
      else
  
        @metadata[item]
      end
    end
    self.metadata_key = items.join('-')
  end

  def create_partial_file()

    partial_file = File.new(self.partial_file_name, 'w')
    partial_file.close()
  end

  def normalize_name(name)

    nm = name.dup
    nm.gsub!('/ ', '_')
    nm.gsub!('(', '')
    nm.gsub!(')', '')
    nm.gsub!(' ', '_')
    nm.gsub!('/', '_')
    nm.gsub!('-', '_')
    return nm
  end
end