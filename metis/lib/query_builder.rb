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
  end
end