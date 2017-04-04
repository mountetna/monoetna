module FileModel

  class File < Sequel::Model

#    def to_hash()
#
#      {
#        :original_name=> original_name,
#        :file_name=> file_name,
#        :file_size=> file_size,
#        :group_name=> group_name,
#        :project_name=> project_name,
#        :start_upload=> start_upload,
#        :finish_upload=> finish_upload,
#        :upload_by=> upload_by,
#        :hashing_algorithm=> hashing_algorithm,
#        :hash=> hash
#      }
#    end
  end

  class Upload < Sequel::Model

    one_to_one :file

#def to_hash()

#  {
#    :original_name=> file.original_name,
#    :file_name=> file.file_name,
#    :file_size=> file.file_size,
#    :group_name=> file.group_name,
#    :project_name=> file.project_name,
#    :start_upload=> file.start_upload,
#    :finish_upload=> file.finish_upload,
#    :upload_by=> file.upload_by,
#    :hashing_algorithm=> file.hashing_algorithm,
#    :hash=> file.hash,

#    :current_byte_position=> current_byte_position,
#    :current_blob_size=> current_blob_size,
#    :next_blob_size=> next_blob_size,
#    :next_blob_hash=> next_blob_hash,
#    :status=> status
#  }
#end
  end
end