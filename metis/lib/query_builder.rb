require 'date'
require_relative './globber'

class Metis
  class QueryBuilder

    def initialize(base_query, params, limit = nil, offset = 0)
      @base_query = base_query
      @params = params
      @limit = limit
      @offset = offset
    end

    def build
      result = @base_query
      @params.each do |param|
        next unless param.key?(:attribute) && param.key?(:predicate) && param.key?(:value)

        # Modify query behavior based on "type" of attribute
        case param[:attribute]
        when 'created_at', 'updated_at'
          result = result.where{|o|
            Sequel[param[:attribute].to_sym].method(param[:predicate]).(Time.iso8601(param[:value]))
          }
        when 'name'
          case param[:predicate]
          when '=~'
            value = param[:value]
            if is_regex(value)
              result = result.where([[Sequel[model_name_attribute], to_regex(value)]])
            else
              result = result.where(Sequel.like(model_name_attribute, value))
            end
          when '='
            result = result.where([[Sequel[model_name_attribute], param[:value]]])
          when 'glob'
            globber = Metis::Globber.new(param[:value], is_file_query)

            id_param = is_file_query ? :folder_id : :id

            result = result.where(
              [[Sequel[id_param], globber.folder_path_ids]]).where(
              Sequel.like(model_name_attribute, globber.sql_search_string))
          end
        when 'folder_id'
          case param[:predicate]
          when 'in'
            raise Metis::QueryError.new("value must be Array of integers for a folder_id search.") unless param[:value].is_a?(Array) && param[:value].all? {|v| v.is_a?(Integer)}
            result = result.where([[:folder_id, param[:value]]])
          end
        end
      end

      result = result.order(:updated_at, :id)
      result = result.offset(@offset)
      result = result.limit(@limit) unless @limit.nil?
      result
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