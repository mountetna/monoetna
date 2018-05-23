
class Metis
  class File < Sequel::Model
    plugin :timestamps, update_on_create: true
    many_to_one :bucket
    many_to_one :folder, class: Metis::File

    FILENAME_MATCH=/[^<>:;,?"*\|\/\x00-\x1f]+/x

    FILEPATH_MATCH=%r!
      \A
        (?<folders>(#{FILENAME_MATCH.source}/)*)

        #{FILENAME_MATCH.source}
      \z
    !x

    def self.valid_file_name?(file_name)
      !!(file_name =~ FILEPATH_MATCH)
    end

    def self.folder_name(file_name)
      FILEPATH_MATCH.match(file_name)[:folders].tap do |folders|
        # no trailing slash in folder_name
        return folders.empty? ? nil : folders.sub(%r!/$!,'')
      end
    end

    def self.safe_file_name(file_name)
      file_name.unpack('H*').first
    end

    def self.md5(path)
      # use md5sum to avoid reading blob
      %x{ md5sum '#{path}' }.split.first
    end

    def self.upload_url(request, project_name, bucket_name, file_name)
      hmac_url(
        'POST',
        request.host,
        Metis::Server.route_path(
          request,
          :upload,
          project_name: project_name,
          bucket_name: bucket_name,
          file_name: file_name
        ),
        Metis.instance.config(:upload_expiration)
      )
    end

    def self.download_url(request, project_name, bucket_name, file_name)
      hmac_url(
        'GET',
        request.host,
        Metis::Server.route_path(
          request,
          :download,
          project_name: project_name,
          bucket_name: bucket_name,
          file_name: file_name
        ),
        Metis.instance.config(:download_expiration)
      )
    end

    def self.author(user)
      [ user.email, "#{user.first} #{user.last}" ].join('|')
    end

    private

    def self.hmac_url(method, host, path, expiration=0)

      URI::HTTPS.build(
        Etna::Hmac.new(
          Metis.instance,
          method: method,
          host: host,
          path: path,
          expiration: (Time.now + expiration).iso8601,
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

    def folder?
      self.is_folder
    end

    def read_only?
      read_only
    end

    def location
      ::File.expand_path(::File.join(
        Metis.instance.project_path(project_name),
        bucket.name,
        Metis::File.safe_file_name(
          file_name
        )
      ))
    end

    def to_hash(request=nil)
      result = {
        file_name: file_name,
        project_name: project_name,
        updated_at: updated_at,
        created_at: created_at,
        author: author,
      }
      if folder?
        return result.merge(to_hash_folder(request))
      else
        return result.merge(to_hash_file(request))
      end
    end

    def to_hash_file(request)
      {
        file_hash: file_hash,
        size: actual_size,
        download_url: request ? Metis::File.download_url(
          request,
          project_name,
          bucket.name,
          file_name
        ) : nil
      }
    end

    def to_hash_folder(request)
      {
        is_folder: true
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
