module Etna
  module Clients
    class Magma
      # Similar to the Attribute.entry method which converts x-form values to inner json entries, but different in
      # that the format is intended to serve etna csv files.  So while the implementations are similar, it is NOT
      # recommended to consolidate the two.  One of them is a backend implementation detail, the other is a stable
      # interface for csvs.
      # Note: Exceptions returned from this method are intended to
      def self.cast_csv_value(string_value, attr, recip)
        # If a file.
        if string_value.empty?
          return "::blank"
        end

        case attr.attribute_type
        when AttributeType::IDENTIFIER
          string_value
        when AttributeType::LINK
          string_value
        when AttributeType::STRING
          string_value
        when AttributeType::MATRIX
          string_value.split(',').map(&:strip!).map(&method(:to_numeric_or_import_error))
        when AttributeType::PARENT
          string_value
        when AttributeType::CHILD
          string_value
        when AttributeType::INTEGER
          to_numeric_or_import_error(string_value, [Integer])
        when AttributeType::TABLE
          raise TypeError.new("Table attributes like #{attr.name} should be specified through the rows of #{attr.link_model_name} linked via the parent attribute.")
        when AttributeType::COLLECTION
          raise TypeError.new("Collection attributes like #{attr.name} should be specified through the rows of #{attr.link_model_name} linked via the parent attribute.")
        when AttributeType::BOOLEAN
          if ['true', 'y', 't', 'yes', 'false', 'n', 'f', 'no'].include?(string_value.downcase)
            Etna::CsvImporter::COLUMN_AS_BOOLEAN.call(string_value)
          else
            raise ArgumentError.new("#{string_value} is not a valid boolean")
          end
        when AttributeType::DATE_TIME
          begin
            DateTime.parse(string_value).iso8601
          rescue Date::Error
            raise ArgumentError.new("#{string_value} is not a valid date time")
          end
        when AttributeType::MATCH
          parse_csv_match(string_value)
        when AttributeType::FLOAT
          to_numeric_or_import_error(string_value, Float)
        when AttributeType::FILE
        when AttributeType::FILE_COLLECTION
        end
      end

      def self.parse_csv_match(string_value)
        comma = string_value.index(',')
        if comma.nil? || comma == 0 || comma == string_value.length
          raise ArgumentError.new("#{string_value} is not a valid match, should be comma separated containing type, value")
        end

        type = string_value.slice(0, comma).strip!
        value = string_value.slice((comma + 1)..-1).strip!

        case type
        when "Array"
          begin
            value = value.split(',').map(&:strip!).map { |v| JSON.parse(v) }
          rescue JSON::ParserError
            raise ArgumentError.new("#{string_value} is not a valid Array match, value after first comma should be encoded as a JSON array")
          end
        when "Regexp"
          begin
            value = Regexp.new(value).source
          rescue RegexpError
            raise ArgumentError.new("#{string_value} is not a valid Regexp match, value after first comma should be a valid regex, but it is not.")
          end
        else
          raise ArgumentError.new("#{string_value} is not a valid match, first value must be either Regexp or Array.")
        end

        {type: type, value: value}
      end

      def self.to_numeric_or_import_error(string_value, numeric_types = [Integer, Float])
        numeric_types.each do |type|
          begin
            return type.call(string_value)
          rescue ArgumentError => e
            next
          end
        end

        raise ArgumentError.new("#{string_value} is not a valid #{numeric_types.map(&:name).join(' or ')}")
      end
    end
  end
end