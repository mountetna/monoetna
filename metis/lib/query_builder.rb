require 'date'

class Metis
  class QueryBuilder

    def initialize(base_query, params)
      @base_query = base_query
      @params = params
    end

    def build
      @params.each do |param|
        next unless param.key?(:attribute) && param.key?(:predicate) && param.key?(:value)

        # Modify query behavior based on "type" of attribute
        case param[:attribute]
        when 'created_at', 'updated_at'
          @base_query = @base_query.where{|o|
            Sequel[param[:attribute].to_sym].method(param[:predicate]).(Time.iso8601(param[:value]))
          }
        when 'name'
          case param[:predicate]
          when '=~'
            value = param[:value]
            if is_regex(value)
              @base_query = @base_query.where([[Sequel[model_name_attribute], to_regex(value)]])
            else
              @base_query = @base_query.where(Sequel.like(model_name_attribute, value))
            end
          when '='
            @base_query = @base_query.where([[Sequel[model_name_attribute], param[:value]]])
          when 'glob'
            glob_parts = param[:value].split('/')

            folder_path_ids = get_folder_path_ids(glob_parts)

            if is_file_query
              id_param = :folder_id
            else
              id_param = :id
            end
            @base_query = @base_query.where(
              [[Sequel[id_param], folder_path_ids]]).where(
              Sequel.like(model_name_attribute, likeify_glob(glob_parts[-1])))
          end
        end
      end
      @base_query
    end

    private

    def to_regex(value)
      return Regexp.new value[1..-2] if value[0] == '/' && value[-1] == '/'

      # otherwise has options
      split_value = value.split('/')
      options = split_value[-1]

      Regexp.new(
        split_value[1],
        (options.include?('i') ? Regexp::IGNORECASE : 0) | (options.include?('e') ? Regexp::MULTILINE : 0) | (options.include?('x') ? Regexp::EXTENDED : 0))
    end

    def is_regex(value)
      # '/foo/' OR '/foo/ix'
      (value[0] == '/' && value[-1] == '/') ||
      (value.split('/').length == 3 && value.split('/')[-1] =~ /[ixe]+/)
    end

    def is_file_query
      @base_query.model == Metis::File
    end

    def model_name_attribute
      is_file_query ? :file_name : :folder_name
    end

    def likeify_glob(glob_string)
      glob_string.gsub('*', '%')
    end

    def recursive_glob(glob_parts)
      (glob_parts.length == 3 && glob_parts[1] == '**') ||
      (glob_parts.length == 2 && glob_parts[1] == '*')
    end

    def depth_one_glob(glob_parts)
      return glob_parts.length == 3 && glob_parts[1] == '*' if is_file_query

      return glob_parts.length == 2 && glob_parts[1].include?('*')
    end

    def get_folder_path_ids(glob_parts)
      # Maybe we only support some subset of pseudo-glob syntax
      folder_name = glob_parts[0]
      if folder_name.include?('*')
        # MVIR1*/
        folders = Metis::Folder.where(Sequel.like(:folder_name, likeify_glob(folder_name))).all
        return folders.map { |f| f.child_folders.map { |f2| f2.id } }.flatten

        # Are there other glob-type folder queries here, like MVIR*/foo* that we should support?
      else
        folder = Metis::Folder.where(folder_name: folder_name).first

        # foo/**/*.txt
        return folder.child_folders.map { |f| f.id } if recursive_glob(glob_parts)

        # 1 level deep glob, like foo/*/*.txt?
        # or foo/bar*
        return folder.child_folders.select { |f|
          f.folder_id == folder.id }.map { |f| f.id } if depth_one_glob(glob_parts)

        # *.txt in the root directory
        return [folder.id]
      end
    end
  end
end