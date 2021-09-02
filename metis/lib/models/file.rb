
class Metis
  class File < Sequel::Model
    plugin :timestamps, update_on_create: true

    many_to_one :bucket
    many_to_one :folder
    many_to_one :data_block

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

    def self.valid_file_name?(file_name)
      !!(file_name =~ /\A#{FILENAME_MATCH.source}\z/)
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

      self.from_folder(bucket, folder, file_name)
    end

    def self.from_folder(bucket, folder, file_name)
      return self.where(
        bucket: bucket,
        folder_id: folder ? folder.id : nil,
        file_name: file_name
      ).first
    end

    def self.safe_file_name(file_name)
      file_name.unpack('H*').first
    end

    def self.unsafe_file_name(encoded_file_name)
      [encoded_file_name].pack('H*')
    end

    def self.md5(path)
      # use md5sum to avoid reading blob
      raise 'Invalid path for md5' unless Metis::File.valid_file_path?(path[1..-1])
      %x{ md5sum "#{path}" }.split.first
    end

    def self.upload_url(request, project_name, bucket_name, file_path, user)
      hmac_url(
        method: 'POST',
        host: request.host,
        path: Metis::Server.route_path(
          request,
          :upload,
          project_name: project_name,
          bucket_name: bucket_name,
          file_path: file_path
        ),
        user: user,
        expiration: Metis.instance.config(:upload_expiration)
      )
    end

    def self.download_url(request, project_name, bucket_name, file_path)
      hmac_url(
        method: 'GET',
        host: request.host,
        path: Metis::Server.route_path(
          request,
          :download,
          project_name: project_name,
          bucket_name: bucket_name,
          file_path: file_path
        ),
        expiration: Metis.instance.config(:download_expiration)
      )
    end

    def self.author(user)
      [ user.email, user.name ].join('|')
    end

    def self.exists?(file_name, bucket, parent_folder)
      Metis::File.where(file_name: file_name, bucket: bucket, folder_id: parent_folder ? parent_folder.id : nil).count > 0
    end

    def self.link_to_block(params)
      # Assumes validations have been done, this just executes
      #   the copy and returns the new copy.
      dest_folder_path, dest_file_name = Metis::File.path_parts(params[:dest_file_path])

      dest_bucket = Metis::Bucket.find(
        project_name: params[:dest_project_name] || params[:project_name],
        name: params[:dest_bucket_name]
      )

      dest_folder = Metis::Folder.from_path(dest_bucket, dest_folder_path).last

      if Metis::File.exists?(dest_file_name, dest_bucket, dest_folder)
        old_dest_file = Metis::File.from_path(dest_bucket, params[:dest_file_path])
        old_dest_file.data_block = params[:source_data_block]
        old_dest_file.save
        return old_dest_file
      else
        return Metis::File.create(
          project_name: dest_bucket.project_name,
          file_name: dest_file_name,
          folder_id: dest_folder&.id,
          bucket: dest_bucket,
          author: Metis::File.author(params[:user]),
          data_block: params[:source_data_block]
        )
      end
    end

    private

    def self.hmac_url(method:, host:, path:, user:  nil, expiration: 0)

      URI::HTTPS.build(
        Etna::Hmac.new(
          Metis.instance,
          method: method,
          host: host,
          path: path,
          expiration: (Time.now + expiration).iso8601,
          nonce: Metis.instance.sign.uid,
          id: :metis,
          headers: user ? {
            email: user.email,
            name: user.name.strip
          } : {}
        ).url_params
      )
    end

    public

    one_to_many :uploads

    def file_hash
      data_block.md5_hash
    end

    def read_only?
      read_only
    end

    def to_hash(request: nil, file_path: nil, with_path: true)
      

      params = {
        folder_id: folder_id,
        file_name: file_name,
        project_name: project_name,
        bucket_name: bucket.name,
        updated_at: updated_at.iso8601,
        created_at: created_at.iso8601,
        author: author,
        file_hash: file_hash,
        archive_id: data_block.archive_id,
        read_only: read_only?,
        size: data_block.actual_size,
      }

      if with_path
        file_path ||= self.file_path
        params[:file_path] = file_path
        params[:download_url] = request ? Metis::File.download_url(
          request,
          project_name,
          bucket.name,
          file_path
        ) : nil
      end

      params
    end

    def file_path
      folder ? ::File.join(folder.folder_path, file_name) : file_name
    end

    def has_data?
      data_block.has_data?
    end

    def can_remove?
      !read_only?
    end

    def remove!
      delete
    end

    def protect!
      update(read_only: true)
    end

    def unprotect!
      update(read_only: false)
    end

    def rename!(new_folder, new_file_name, user=nil)
      new_params = {
        file_name: new_file_name,
        folder_id: new_folder ? new_folder.id : nil
      }

      new_params[:author] = Metis::File.author(user) if user

      update(**new_params)
    end

    def update_bucket!(new_bucket)
      raise 'Bucket does not match folder bucket' if folder != nil && folder.bucket_id != new_bucket.id
      update(bucket: new_bucket)
      refresh
    end

    def update_bucket_and_rename!(folder, new_file_name, new_bucket, user=nil)
      raise 'Bucket does not match folder bucket' if folder != nil && folder.bucket_id != new_bucket.id
      new_params = {
        folder_id: folder ? folder.id : nil,
        file_name: new_file_name,
        bucket: new_bucket
      }

      new_params[:author] = Metis::File.author(user) if user

      update(**new_params)
      refresh
    end
  end
end
