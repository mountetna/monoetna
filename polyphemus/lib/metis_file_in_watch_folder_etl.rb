# Returns a set of Metis files, if they are in a watch folder
require_relative "./metis_file_etl"

class Polyphemus
  class MetisFileInWatchFolderCursor < MetisFileEtlCursor
    def initialize(job_name:, project_name:, bucket_name:)
      raise "project_name cannot be nil" if project_name.nil?
      raise "bucket_name cannot be nil" if bucket_name.nil?
      super("#{job_name}_metis_files_in_watch_folder_#{project_name}_#{bucket_name}")
      self[:project_name] = project_name
      self[:bucket_name] = bucket_name
    end

    def reset!
      super { self[:seen_ids] = [] }
    end
  end

  class MetisFileInWatchFolderEtl < MetisFileEtl
    def execute_request(find_request)
      # Only pass along files that are in a watch folder
      files = super(find_request)

      files.select do |file|
        watch_folder_paths(
          find_request.project_name, find_request.bucket_name
        ).include?(::File.dirname(file.file_path))
      end
    end

    def watch_folder_paths(project_name, bucket_name)
      @watch_folder_paths ||= Polyphemus::WatchFolder.where(
        project_name: project_name,
        bucket_name: bucket_name,
      ).all.map do |folder|
        folder.folder_path
      end
    end
  end
end
