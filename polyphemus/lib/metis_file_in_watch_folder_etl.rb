# Returns a set of Metis files, if they are in a watch folder
require_relative "./metis_file_etl"

class Polyphemus
  class MetisFileInWatchFolderEtl < MetisFileEtl
    attr_reader :bucket_watch_configs

    def initialize(metis_client: nil, bucket_watch_configs: [], limit: 20)
      pairs = bucket_watch_configs.map(&:cursor_pair)
      raise "bucket_watch_configs have non unique project/bucket combinations" unless pairs.uniq.length == pairs.length

      @bucket_watch_configs = bucket_watch_configs
      super(
        project_bucket_pairs: pairs,
        metis_client: metis_client,
        limit: limit,
        file_name_params: {},
      )
    end

    def group_by_watch_type(cursor, files)
      watch_folders = watch_folders(cursor)

      {}.tap do |result|
        files.each do |file|
          matching_watches = watch_folders.select { |wf| wf.metis_id == file.folder_id }
          matching_watches.each do |watch|
            matched_files = (result[watch.watch_type] ||= [])
            matched_files << file
          end
        end
      end
    end

    def process(cursor, files)
      group_by_watch_type(cursor, flatten_grouped_files(files)).each do |watch_type, files|
        process_files(cursor, files, watch_type)
      end
    end

    def process_files(cursor, files, watch_type)
      raise "Subclasses should implement process_files"
    end

    def prepare_find_request(cursor, find_request)
      super

      # Here we update the find request with the folder_ids from
      #   the watch folders
      find_request.add_param(Etna::Clients::Metis::FindParam.new(
        type: "file",
        attribute: "folder_id",
        predicate: "in",
        value: watch_folders(cursor).map(&:metis_id),
      ))
    end

    def execute_request(find_request, i)
      raw_files = super

      raw_files.group_by(&:folder_id).map do |folder_id, files|
        MetisFilesInFolder.new(folder_id, files)
      end
    end

    def matching_bucket_config(cursor)
      @bucket_watch_configs.find do |c|
        c.project_name == cursor[:project_name] \
          && c.bucket_name == cursor[:bucket_name]
      end.tap do |bucket_config|
        if bucket_config.nil?
          raise "Could not find a configured bucket config for the cursor, possibly reset the cursor or fix etl config"
        end
      end
    end

    def watch_folders(cursor)
      bucket_config = matching_bucket_config(cursor)

      Polyphemus::WatchFolder.where(
        project_name: cursor[:project_name],
        bucket_name: cursor[:bucket_name],
        watch_type: bucket_config.all_watch_types,
      ).all
    end

    def flatten_grouped_files(grouped_metis_files)
      grouped_metis_files.map(&:files).flatten
    end
  end

  class MetisFilesInFolder
    attr_reader :identifier, :files

    def initialize(identifier, metis_files)
      @identifier = identifier
      @files = metis_files
    end

    def file_paths_hashes
      files.map do |file|
        [file.file_path, file.file_hash]
      end
    end

    def file_path
      files.map(&:file_path).join(",")
    end

    def updated_at
      files.map do |file|
        file.updated_at
      end.minmax.last
    end

    def length
      files.length
    end
  end
end
