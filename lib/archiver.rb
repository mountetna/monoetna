require 'fog/aws'

class Metis
  class Archiver
    def initialize(config)
      @config = config
    end

    def archive(file)
      raise ArgumentError, "No vault defined!" unless @config

      backup = Metis::Backup.where(md5_hash: file.file_hash).first

      if !backup
        description = "md5:#{file.file_hash} #{file.file_path}"

        archive_id = nil
        if file.actual_size == 0
          archive_id = 'zero-byte-file'
          description = 'A zero-byte file'
        else
          archive = vault.archives.create(
            body: ::File.open(file.location),
            multipart_chunk_size: 64*1024*1024,
            description: "md5:#{file.file_hash} #{file.file_path}"
          )

          return nil unless archive

          archive_id = archive.id
        end

        backup = Metis::Backup.create(
          md5_hash: file.file_hash,
          description: description,
          archive_id: archive_id
        )
        backup.save
      end

      file.update(backup_id: backup.id)
    end

    private

    def vault
      glacier.vaults.find { |v|
        v.id == @config[:directory]
      } || glacier.vaults.create(
        id: @config[:directory]
      )
    end

    def glacier
      @glacier ||= Fog::AWS::Glacier.new(
        @config[:credentials]
      )
    end
  end
end
