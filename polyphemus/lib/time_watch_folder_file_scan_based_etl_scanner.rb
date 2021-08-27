# Returns a set of Metis files, if they are in a watch folder
require_relative "./time_scan_based_etl_scanner"

class Polyphemus
  class TimeWatchFolderFileScanBasedEtlScanner < TimeScanBasedEtlScanner
    def find_batch(cursor)
      require "pry"
      binding.pry
      files = super
      files.select do |file|
        watch_folder_paths(
          cursor[:project_name], cursor[:bucket_name]
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
