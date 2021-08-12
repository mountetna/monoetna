class Metis
  class Globber
    def initialize(glob_string, is_file_glob)
      @glob_parts = glob_string.split("/")
      @is_file_glob = is_file_glob
    end

    def sql_search_string
      likeify_glob(@glob_parts[-1])
    end

    def folder_path_ids
      # Maybe we only support some subset of pseudo-glob syntax
      folder_name = @glob_parts[0]
      if folder_name.include?("*")
        # MVIR1*/
        folders = Metis::Folder.where(Sequel.like(:folder_name, likeify_glob(folder_name))).all
        return folders.map { |f| f.id } + folders.map { |f| f.child_folders.map { |f2| f2.id } }.flatten
      else
        folders = Metis::Folder.where(folder_name: folder_name).all

        # foo/**/*.txt
        return folders.map { |f| f.id }.concat(folders.map { |f| f.child_folders.map { |f2| f2.id } }).flatten if recursive_glob

        # 1 level deep glob. For files, like foo/*/*.txt
        # or foo/bar* for directories
        return folders.map { |f|
                 f.child_folders.select { |f2|
                   f2.folder_id == f.id
                 }.map { |f3| f3.id }
               }.flatten if depth_one_glob

        # *.txt in the root directory for files, or foo/*
        return folders.map { |f| f.id }
      end
    end

    private

    def likeify_glob(glob_string)
      glob_string.gsub("*", "%")
    end

    def recursive_glob
      return @glob_parts.length == 3 && @glob_parts[1] == "**" if @is_file_glob

      (@glob_parts.length == 2 && @glob_parts[1] == "*") ||
      (@glob_parts.length == 3 && @glob_parts[1] == "**")
    end

    def depth_one_glob
      return @glob_parts.length == 3 && @glob_parts[1] == "*" if @is_file_glob

      @glob_parts.length == 2 && @glob_parts[1].include?("*")
    end
  end
end
