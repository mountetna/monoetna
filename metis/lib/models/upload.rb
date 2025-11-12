class Metis
  class Upload < Sequel::Model
    many_to_one :bucket

    def to_hash
      [:project_name, :file_name, :author, :current_byte_position, :next_blob_size, :next_blob_hash].map do |s|
        [s, send(s)]
      end.to_h
    end

    def to_json
      to_hash.to_json
    end

    def partial_location
      ::File.expand_path(
          ::File.join(
              Metis.instance.config(:data_path),
              'uploads',
              Digest::MD5.hexdigest("#{metis_uid}-#{id.to_s}")
          )
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
      self.delete_partial!
      delete
    end

    def delete_partial!
      if ::File.exists?(partial_location)
        ::File.delete(partial_location)
      end
    end

    def finish!
      folder_path, new_file_name = Metis::File.path_parts(file_name)

      folder = Metis::Folder.from_path(bucket, folder_path).last

      data_block = Metis::DataBlock.create_from(file_name, partial_location)

      file = Metis::File.find_or_create(
          project_name: project_name,
          file_name: new_file_name,
          folder_id: folder&.id,
          bucket: bucket
      ) do |f|
        f.author = author
        f.data_block = data_block
      end
      file.update(folder: folder, author: author, data_block: data_block)
      Metis::DataBlockLedger.log_create(
        file,
        data_block,
        author
      )
      Metis::DataBlockLedger.log_link(file, data_block, author)
      return file
    end

    def complete?
      file_size == ::File.size(partial_location)
    end

    def self.fetch(params)
      Metis::Upload.find_or_create(
          file_name: params[:file_name],
          bucket: params[:bucket],
          metis_uid: params[:metis_uid],
          project_name: params[:project_name]
      ) do |upload|
        upload.set(
            author: Metis::File.author(params[:user]),
            file_size: 0,
            current_byte_position: 0,
            next_blob_size: -1,
            next_blob_hash: '',
        )
      end
    end

  end
end
