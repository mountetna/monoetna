class Metis
  class Upload < Sequel::Model
    many_to_one :file

    def to_hash
      [ :current_byte_position, :next_blob_size, :next_blob_hash ].map do |s|
        [ s, send(s) ]
      end.to_h.merge(
        project_name: file.project_name,
        file_name: file.file_name
      )
    end

    def to_json
      to_hash.to_json
    end

    def partial_location
      ::File.expand_path(::File.join(
        Metis.instance.project_path(file.project_name),
        'uploads',
        "#{metis_uid}-#{file.file_name}"
      ))
    end

    def append_blob(blob_path)
      # use cat to avoid reading file
      %x{ cat #{blob_path} >> "#{partial_location}" }

      self.update(
        current_byte_position: ::File.size(partial_location)
      )
    end

    def finish!
      # Update the postgres record.
      file.set_file_data(partial_location)
    end

    def complete?
      file_size == ::File.size(partial_location)
    end

    def blob_valid?(next_blob_path)
      # next_blob_hash and _size are the expected
      # content hash and size of the blob

      return (
        Metis::File.md5(next_blob_path) == next_blob_hash &&
        ::File.size(next_blob_path) == next_blob_size
      )
    end
  end
end
