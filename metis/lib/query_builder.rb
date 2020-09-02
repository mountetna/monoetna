require 'date'

class Metis
  class QueryBuilder

    attr_reader :base_query

    def initialize(base_query, params)
      @base_query = base_query
      @params = params
    end

    def build
      @params.each do |param|
        next unless param.key?(:attribute) && param.key?(:predicate) && param.key?(:value)

        # Modify query behavior based on "type" of attribute....but no way
        #   to know that from the user params?
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

    def model_name_attribute
      @base_query.model == Metis::File ? :file_name : :folder_name
    end
  end
end