require_relative "../tsv_writer_base"
require_relative "./question_format"

class Magma
  class QueryTSVWriter < Magma::TSVWriterBase
    def initialize(project_name, question, opts = {})
      @project_name = project_name
      @question = question

      validate_columns(opts[:user_columns]) if opts[:user_columns]

      @user_columns = opts[:user_columns]

      super(opts)
    end

    def write_tsv(&block)
      @transpose ?
        transpose_records(&block) :
        standard_records(&block)
    end

    def standard_records
      yield tsv_header

      @question.each_page_answer do |answer|
        yield to_tsv(answer)
      end
    end

    def transpose_records
      # Collect all data, then yield the transposed results
      [].tap do |buffer|
        buffer << tsv_row_to_array(tsv_header)

        buffer.concat(tsv_records_to_array(to_tsv(@question.answer)))
      end.transpose.each do |row|
        yield array_to_tsv_row(row)
      end
    end

    private

    def validate_columns(user_columns)
      raise TSVError.new("user_columns array must be #{model_attr_headers.length} elements long") unless model_attr_headers.length == user_columns.length
    end

    def path_to_value(search_array, target_column, current_path: [], starting_index: nil)
      # Given search_array as a nested set of arrays,
      #   finds the index-path to get to a specific value, starting with
      #   the starting_index if provided.
      # i.e. search_array is [labors::labor#name, [labors::labor#number, labors::monster#name]]
      #   the path to labors::monster#name is [1, 1]
      return [] unless search_array
      search_array = [search_array] unless search_array.is_a?(Array)

      return current_path.concat([starting_index]) if starting_index && search_array[starting_index] == target_column.header

      direct_index = search_array.find_index(target_column.header)
      return current_path.concat([direct_index]) unless direct_index.nil?

      refined_search_array = starting_index.nil? ? search_array : [search_array[starting_index]]

      refined_search_array.each.with_index do |element, index|
        if element.is_a?(Array)
          temp_path = path_to_value(
            element,
            target_column,
            current_path: current_path.concat([starting_index.nil? ? index : starting_index]),
          )

          return temp_path unless temp_path.empty?
        end
      end

      []
    end

    def model_attr_headers
      # "raw" headers that reference only the model + attribute names
      @model_attr_headers ||= @question.columns.map do |col|
        TSVHeader.new(@project_name, col)
      end
    end

    def rename(header, index)
      @user_columns ? @user_columns[index] : header.header
    end

    def matrix_columns(header, index)
      path_to_matrix_attribute = path_to_value(
        @question.format[1],
        header,
        starting_index: index - 1,
      )

      @question.format[1].dig(*(path_to_matrix_attribute.slice(0..-2).concat([1])))
    end

    def expand(header, index)
      renamed_header = rename(header, index)
      matrix_columns(header, index).map do |col|
        "#{renamed_header}.#{col}"
      end
    end

    def tsv_header
      # Start with the raw, internal headers.
      # If the user supplies a :user_columns option, in
      #   which case, rename according to the :display_label
      # Expand matrix headers if necessary
      model_attr_headers.map.with_index do |model_attr_header, index|
        @expand_matrices && model_attr_header.matrix? ?
          expand(model_attr_header, index) :
          rename(model_attr_header, index)
      end.flatten.join("\t") + "\n"
    end

    def to_tsv(records)
      CSV.generate(col_sep: "\t") do |csv|
        records.map do |record|
          csv << [].tap do |row|
            model_attr_headers.each.with_index do |tsv_column, index|
              column_index = index - 1

              if index == 0
                row << record.identifier
                next
              elsif non_nested_single_model_query?
                # In this simple use case, we just grab the entire
                #   answer portion
                row << extract_data(record.data)
                next
              else
                attribute_data = record.data[column_index]

                if attribute_data.data.is_a?(Magma::MatrixPredicate::MatrixValue)
                  value = JSON.parse(attribute_data.data.to_json)
                elsif attribute_data.is_a?(Magma::SimpleAnswerBase)
                  value = extract_data(attribute_data.data)
                else
                  value = attribute_data.aggregated_values(tsv_column.array?)
                end
              end

              if expand_matrix_data(tsv_column.matrix?)
                row = row.concat(value.nil? ?
                  Array.new(matrix_columns(tsv_column, column_index).length) { nil } :
                  value)
              else
                row << value
              end
            end
          end
        end
      end
    end

    def non_nested_single_model_query?
      # Simple edge case when query is something like
      #   [model, ::all, attribute]
      @question.format.length == 2 &&
      @question.format.last.is_a?(String)
    end

    def expand_matrix_data(is_matrix_column)
      @expand_matrices && is_matrix_column
    end
  end

  class TSVHeader < Magma::QuestionColumnBase
    attr_reader :header

    def initialize(project_name, header)
      @header = header
      @project_name = project_name
    end

    def matrix?
      is_matrix?(@header)
    end

    def array?
      is_collection?(@header)
    end
  end
end
