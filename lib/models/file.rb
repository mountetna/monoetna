class Metis
  class File < Sequel::Model
    def self.md5(path)
      # use md5sum to avoid reading blob
      %x{ md5sum '#{path}' }.split.first
    end

    def self.upload_url(request, project_name, file_name)
      hmac_url(
        'POST',
        request.host,
        Metis::Server.route_path(
          request,
          :upload,
          project_name: project_name,
          file_name: file_name
        ),
        (Time.now + Metis.instance.config(:upload_expiration)).iso8601
      )
    end

    def self.download_url(request, project_name, file_name)
      hmac_url(
        'GET',
        request.host,
        Metis::Server.route_path(
          request,
          :download,
          project_name: project_name,
          file_name: file_name
        ),
        (Time.now + Metis.instance.config(:download_expiration)).iso8601
      )
    end

    private

    def self.hmac_url(method, host, path, expiration)

      URI::HTTPS.build(
        Etna::Hmac.new(
          Metis.instance,
          method: method,
          host: host,
          path: path,
          expiration: expiration,
          nonce: Metis.instance.sign.uid,
          id: :metis,
          headers: { }
        ).url_params
      )
    end

    public

    one_to_many :uploads

    def self.has_file?(project_name, file_name)
      file = self.where(project_name: project_name, file_name: file_name).first

      return file && file.has_data?
    end

    def compute_hash!
      if has_data?
        self.file_hash = Metis::File.md5(location)
        save
      end
    end

    def has_data?
      file_name && ::File.exists?(location)
    end

    def location
      ::File.expand_path(::File.join(
        Metis.instance.project_path(project_name), file_name
      ))
    end

    def to_hash(request=nil)
      {
        file_name: file_name,
        project_name: project_name,
        original_name: original_name,
        size: actual_size,
        file_hash: file_hash,
        download_url: request ? Metis::File.download_url(
          request,
          project_name,
          file_name
        ) : nil
      }
    end

    def actual_size
      has_data? ? ::File.size(location) : nil
    end

    def set_file_data(file_path)
      # Rename the existing file.
      ::File.rename(
        file_path,
        location
      )

      # update the hash
      self.update(file_hash: Metis::File.md5(location))
    end
  end
end
