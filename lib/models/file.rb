
class Metis
  class File < Sequel::Model
    plugin :timestamps, update_on_create: true
    plugin :rcte_tree,
      parent: {name: :folder},
      ancestors: {name: :folder_path}

    many_to_one :bucket
    many_to_one :folder

    FILENAME_MATCH=/[^<>:;,?"*\|\/\x00-\x1f]+/x

    FILEPATH_MATCH=%r!
      \A
        (?<folders>(#{FILENAME_MATCH.source}/)*)

        (?<name>#{FILENAME_MATCH.source})
      \z
    !x

    def self.valid_file_path?(file_path)
      !!(file_path =~ FILEPATH_MATCH)
    end

    def self.path_parts(file_path)
      FILEPATH_MATCH.match(file_path) do |match|
        return [
          # folder part, with no trailing slash
          match[:folders].empty? ? nil : match[:folders].sub(%r!/$!,''),

          # name part, with no trailing slash
          match[:name].sub(%r!/$!,'')
        ]
      end
    end

    def self.from_path(bucket, file_path)
      folder_path, file_name = path_parts(file_path)

      folder = folder_path && Metis::Folder.from_path(bucket, folder_path).last

      return self.where(
        bucket: bucket,
        folder_id: folder ? folder.id : nil,
        file_name: file_name
      ).first
    end

    def self.safe_file_name(file_name)
      file_name.unpack('H*').first
    end

    def self.md5(path)
      # use md5sum to avoid reading blob
      %x{ md5sum '#{path}' }.split.first
    end

    def self.upload_url(request, project_name, bucket_name, file_path)
      hmac_url(
        'POST',
        request.host,
        Metis::Server.route_path(
          request,
          :upload,
          project_name: project_name,
          bucket_name: bucket_name,
          file_path: file_path
        ),
        Metis.instance.config(:upload_expiration)
      )
    end

    def self.download_url(request, project_name, bucket_name, file_path)
      hmac_url(
        'GET',
        request.host,
        Metis::Server.route_path(
          request,
          :download,
          project_name: project_name,
          bucket_name: bucket_name,
          file_path: file_path
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

    def read_only?
      read_only
    end

    def location
      ::File.expand_path(
        ::File.join(
          folder ? folder.location : bucket.location,
          Metis::File.safe_file_name(file_name)
        )
      )
    end

    def to_hash(request=nil)
      result = {
        file_name: file_name,
        project_name: project_name,
        bucket_name: bucket.name,
        updated_at: updated_at,
        created_at: created_at,
        author: author,
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
