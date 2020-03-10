require 'fog/aws'

class Metis
  class Archiver
    def initialize(config)
      @config = config
    end

    def archive(data_block)
      raise ArgumentError, "No vault defined!" unless @config

      archive_id = nil

      if data_block.actual_size == 0
        archive_id = 'zero-byte-file'
      else
        archive = vault.archives.create(
          body: ::File.open(data_block.location),
          multipart_chunk_size: 64*1024*1024,
          description: "md5:#{data_block.md5_hash} #{data_block.description}"
        )

        return nil unless archive

        archive_id = archive.id
      end

      data_block.archive_id = archive_id
      data_block.save
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
