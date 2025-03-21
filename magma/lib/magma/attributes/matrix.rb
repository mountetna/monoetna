require 'set'
require 'json'

class Magma
  class MatrixJsonError < StandardError
  end
  class MatrixAttribute < Attribute
    def database_type
      :json
    end

    def entry(value, loader)
      [ column_name, value.to_json ]
    end

    class Validation < Magma::Validation::Attribute::BaseAttributeValidation
      def validate(value, &block)
        # nil is a valid value
        return if value.nil?

        # it must be an array of numbers
        yield "Matrix value is not an array of numbers" unless value.is_a?(Array) && value.all?{|v| v.is_a?(Numeric)}
        yield "Improper matrix row size" unless validation_object.options.size == value.size
      end
    end

    def revision_to_loader(record_name, new_value)
      [ name, new_value ]
    end

    def revision_to_payload(record_name, new_value, loader)
      [ name, new_value ]
    end

    def expand(matrix_value)
      # Provide an embedded data frame inside of the TSV
      #   that includes the selected columns names + values.

      # matrix_value.to_json is actually a JSON string.
      JSON.parse(matrix_value.to_json)
    end

    def query_to_tsv(value)
      value.to_json
    end

    def matrix_row_json(data, column_names)
      # may not have identifier of the matrix model, due to
      #   subselects, so we don't cache the data right now.
      return column_names ? Array.new(column_names.length, nil).to_json : null_row_json unless data

      if column_names
        indexes = column_indexes(column_names)
        data ? data.values_at(
          *indexes
        ).to_json : indexes.map{nil}.to_json
      else
        data.to_json
      end
    end

    private

    def null_row_json
      @null_row_json ||= validation_object.options.map{nil}.to_json
    end

    def column_indexes(names)
      @column_indexes ||= validation_object.options.map.with_index{|name,i| [ name, i ]}.to_h

      @column_indexes.values_at(*names)
    end
  end
end
