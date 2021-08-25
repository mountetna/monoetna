class Polyphemus
  class WatchFolder < Sequel::Model
    plugin :timestamps, update_on_create: true

    def to_hash
      {
        project_name: project_name,
        bucket_name: bucket_name,
        folder_path: folder_path,
      }
    end

    def full_path
      "/#{project_name}/#{bucket_name}/#{folder_path}"
    end
  end
end
