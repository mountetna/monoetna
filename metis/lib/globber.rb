class Metis
  class Globber
    STAR_MAGIC=SecureRandom.hex
    BRACE1_MAGIC=SecureRandom.hex
    BRACE2_MAGIC=SecureRandom.hex

    def initialize(bucket, glob_string, is_file_glob)
      @bucket = bucket
      @glob_parts = glob_string.split("/")
      @is_file_glob = is_file_glob
    end

    def sql_search_string
      likeify_glob(@glob_parts[-1])
    end

    def folder_path_ids
      @bucket.folder_tree.paths.select do |id, path|
        folder_match =~ path
      end.keys
    end

    private

    def folder_glob
      "**/" + @glob_parts[0...-1].join('/')
    end

    def folder_match
      match = folder_glob

      # block out glob syntax
      match.gsub!(%r!\*!, STAR_MAGIC)
      match.gsub!(%r!\{!, BRACE1_MAGIC)
      match.gsub!(%r!\}!, BRACE2_MAGIC)

      # escape remaining regexp symbols
      match = Regexp.escape(match)

      # **/ matches any number of directories
      match.gsub!( %r`(?<=/|^)#{STAR_MAGIC*2}/`, '(?:([^/]*(/|$))*)')

      # /** matches any number of directories
      match.gsub!( %r`/#{STAR_MAGIC*2}$`, '(?:(/[^/]*)*)')

      # * matches any non-slash including nothing
      match.gsub!( %r!#{STAR_MAGIC}!, '[^/]*')

      # {a,b,c,d} matches a group
      match.gsub!(%r!#{BRACE1_MAGIC}((?:[^,]+,)+[^,]+)#{BRACE2_MAGIC}!) do |m|
        "(#{Regexp.last_match[1].sub(',','|')})"
      end

      # terminate - the remaining match comes from sql_search_string
      match = "^#{match}$"

      Regexp.new(match)
    end

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
