class Metis
  class FolderTree
    def initialize(bucket)
      @bucket = bucket
    end

    def paths
      @paths ||= add_folder_paths(nil, nil, []).to_h
    end

    private

    def add_folder_paths(id, root, paths)
      folders[id].each do |folder|
        folder_path = ::File.join(*[root, folder[:folder_name]].compact)
        paths.push([ folder[:id], folder_path ])
        if folders[folder[:id]]
          add_folder_paths(
            folder[:id],
            folder_path,
            paths
          )
        end
      end

      paths
    end
 
    def folders
      @folders ||= begin
        heading = [:folder_name, :id, :folder_id]
        Metis::Folder.where(
          bucket: @bucket
        ).select_map(heading).map do |row|
            heading.zip(row).to_h
        end.group_by do |folder|
          folder[:folder_id]
        end
      end
    end
  end
end
