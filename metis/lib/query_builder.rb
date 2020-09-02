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
            @base_query = @base_query.where(Sequel.like(model_name_attribute, param[:value]))
          end
        end
      end
      @base_query
    end

    private

    def is_file_query
      @base_query.model == Metis::File
    end

    def model_name_attribute
      is_file_query ? :file_name : :folder_name
    end
  end
end