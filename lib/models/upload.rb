class Metis
  class Upload < Sequel::Model
    many_to_one :file

    def to_json
      [ :status, :current_byte_position, :current_blob_size, :next_blob_size, :next_blob_hash ].map do |s|
        [ s, send(s) ]
      end.to_h.merge(
        project_name: file.project_name,
        file_name: file.file_name
      ).to_json
    end
  end
end
