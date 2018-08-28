require 'fog/aws'

class Metis
  class Backup
    def initialize(config)
      @config = config
      @glacier = {}
    end

    def archive(project_name, file)
      project_name = project_name.to_sym
      raise ArgumentError, "No vault defined for project #{project_name}!" unless @config.has_key?(project_name)

      archive = vault(project_name).archives.create(
        body: ::File.open(file.location),
        multipart_chunk_size: 64*1024*1024,
        description: "md5:#{file.file_hash} #{file.file_path}"
      )

      file.update(archive_id: archive.id) if archive
    end

    private

    def vault(project_name)
      glacier(project_name).vaults.find { |v|
        v.id == @config[project_name][:directory]
      } || glacier(project_name).vaults.create(
        id: @config[project_name][:directory]
      )
    end

    def glacier(project_name)
      @glacier[project_name] ||= Fog::AWS::Glacier.new(
        @config[project_name][:credentials]
      )
    end
  end
end
