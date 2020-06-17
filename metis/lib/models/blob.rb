class Metis
  class Blob
    def initialize(blob_data)
      @tempfile = blob_data[:tempfile]

      raise 'Invalid blob!' unless @tempfile.is_a?(Tempfile)
    end

    def continues?(upload)
      hash == upload.next_blob_hash && size == upload.next_blob_size
    end

    def path
      @tempfile.path
    end

    private

    def size
      ::File.size(path)
    end

    def hash
      Metis::File.md5(path)
    end
  end
end
