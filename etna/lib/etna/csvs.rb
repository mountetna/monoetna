module Etna
  class CsvImporter
    def initialize(
        strip: true,
        filter_empties: true,
        &row_formatter
    )
      # Removes any columns from a row that are empty strings.  Allows for some simpler 'empty' processing.  This occurs
      # post stripping as well.
      @filter_empties = filter_empties
      @row_formatter = row_formatter
    end

    COLUMN_AS_BOOLEAN = -> (s) { ['true', 't', 'y', 'yes'].include?(s&.downcase) }

    def feed_rows_into(enumerator, filename: nil, input_io: nil, &block)
      begin
        enumerator.next
      rescue StopIteration
        return
      end

      context = {}
      each_csv_row(filename: filename, input_io: input_io) do |row, lineno|
        p row
        block.call(row, lineno) if block_given?
        enumerator.feed(NestedRowProcessor.new(row, lineno, context))

        begin
          enumerator.next
        rescue StopIteration
          return
        end
      end
    end

    def each_csv_row(filename: nil, input_io: nil, &block)
      if input_io.nil?
        unless filename.nil?
          File.open(filename, 'r') do |io|
            return each_csv_row(input_io: io, &block)
          end
        end
      end

      lineno = 1
      header_mappings = nil
      CSV.parse(input_io) do |row|
        if header_mappings.nil?
          header_mappings = row.each_with_index.map { |heading, index| [index, heading&.to_sym] }.to_h
          next
        end

        lineno += 1

        row_as_hash = {}
        row.each_with_index do |value, index|
          if header_mappings.include?(index) && !(header = header_mappings[index]).nil?
            row_as_hash[header] = value
          else
            row_as_hash[index] = value
          end
        end

        row = row_as_hash

        row.keys.each { |k| row[k].strip! if row[k] =~ /^\s+$/ } if @strip
        row.select! { |k, v| !v.empty? } if @filter_empties
        @row_formatter.call(row) unless @row_formatter.nil?
        yield row, lineno if block_given?
      end
    end

    def replace_row_column(row, column, &block)
      if !row[column].nil? || block.arity == 0
        row[column] = yield row[column]
      end
    end

    class ImportError < StandardError
      attr_accessor :lineno

      def initialize(msg, lineno = nil)
        @lineno = lineno
        super(msg)
      end

      def message
        "line #{lineno}: #{super}"
      end
    end

    class NestedRowProcessor
      attr_reader :row, :lineno, :context, :errors

      def initialize(row, lineno, context)
        @row = row
        @lineno = lineno
        @context = context
        @errors = []

        # If a parent context changes, all child contexts are invalidated.  But since parent contexts are changed
        # before the relationship of child contexts are declared, we have to track that so that when a child context
        # dependency is declared we can clear it based on wether dependencies have changed.
        @changed = {}
      end

      def watch(column, *dependencies, &block)
        prev = @context[column]
        process(column, *dependencies) do |*args|
          block.call(prev, *args)
        end
      end

      def process(column, *dependencies, &block)
        if dependencies.any? { |p| @changed.include?(p) }
          @changed[column] = true
          @context[column] = nil
        end

        return false if !row.include?(column) || (next_val = row[column]).nil?
        @changed[column] = true

        context_values = dependencies.map do |p|
          if @context[p].nil?
            raise ImportError.new("Found a #{column} value, but no previous #{p} had been given!", @lineno)
          end

          @context[p]
        end

        begin
          if block_given?
            next_val = yield next_val, *context_values
          end
        rescue ImportError => e
          e.lineno = @lineno
          raise e
        end

        @context[column] = next_val
        true
      end
    end
  end

  class CsvExporter
    # column_headers should be an array of symbols, mapping the column heading names and ordering to export
    # column_serializer is an optional block that takes column (string), column_value (string) and should
    #   return a string representation of column_value to write to the csv.  By default, when nil, the exporter
    #   will attempt to convert the value to a string via to_s or simply write an empty string for nil.

    attr_reader :column_headers
    def initialize(column_headers, &column_serializer)
      @column_headers = column_headers
      @column_serializer = column_serializer
    end

    def header_row
      @column_headers.map(&:to_s)
    end

    def map_column_value(column, column_value)
      @column_serializer&.call(column, column_value) || column_value&.to_s || ''
    end

    def row_from_columns(columns)
      @column_headers.map { |c| self.map_column_value(c, columns[c] || '') }
    end

    def with_row_writeable(filename: nil, output_io: nil, &block)
      if output_io.nil? && !filename.nil?
        File.open(filename, 'w') do |io|
          return with_row_writeable(output_io: io, &block)
        end
      end

      writeable = self.class::RowWriteable.new(self, CSV.new(output_io))
      yield writeable
      writeable.ensure_headers
      nil
    end

    class RowWriteable
      def initialize(exporter, csv)
        @exporter = exporter
        @csv = csv
        @written_headers = false
      end

      def <<(columns)
        self.write_row(columns, nil)
      end

      def write(**columns)
        self.write_row(columns, nil)
      end

      def write_row(columns, row = nil)
        columns_row = @exporter.row_from_columns(columns)
        columns_row += row if row
        self.ensure_headers
        @csv <<  columns_row
      end

      def ensure_headers
        @csv << @exporter.header_row unless @written_headers
        @written_headers = true
      end
    end
  end
end