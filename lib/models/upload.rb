class Metis
  class Upload < Sequel::Model
    many_to_one :bucket

    def to_hash
      [ :project_name, :file_name, :author, :current_byte_position, :next_blob_size, :next_blob_hash ].map do |s|
        [ s, send(s) ]
      end.to_h
    end

    def to_json
      to_hash.to_json
    end

    def self.ensure_upload_dir!(project_name)
      location = upload_location(project_name)
      ::FileUtils.mkdir_p(location) unless ::File.exists?(location)
    end

    def self.upload_location(project_name, file_name=nil)
      ::File.expand_path(::File.join(
        *[
          Metis.instance.project_path(project_name),
          'uploads',
          file_name
        ].compact
      ))
    end

    def partial_location
      Metis::Upload.upload_location(
        project_name,
        Digest::MD5.hexdigest("#{metis_uid}-#{id.to_s}")
      )
    end

    def append_blob(blob, next_blob_size, next_blob_hash)

      # use cat to avoid reading file
      # Verify that partial_location is well-formed
      raise 'Appending to invalid location' unless Metis::File.valid_file_path?(partial_location[1..-1])
      %x{ cat #{blob.path} >> "#{partial_location}" }

      self.update(
        current_byte_position: ::File.size(partial_location),
        next_blob_size: next_blob_size,
        next_blob_hash: next_blob_hash
      )
    end

    def delete_with_partial!
      if ::File.exists?(partial_location)
        ::File.delete(partial_location)
      end
      delete
    end

    def finish!
      folder_path, new_file_name = Metis::File.path_parts(file_name)

      folder = Metis::Folder.from_path(bucket, folder_path).last

      file = Metis::File.find_or_create(
        project_name: project_name,
        file_name: new_file_name,
        folder_id: folder ? folder.id : nil,
        bucket: bucket
      ) do |f|
        f.author = author
      end

      file.update(folder: folder, author: author)

      file.set_file_data(partial_location)

      return file
    end

    def complete?
      file_size == ::File.size(partial_location)
    end
  end
end
