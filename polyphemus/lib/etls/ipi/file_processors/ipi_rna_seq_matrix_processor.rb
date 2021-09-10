require_relative "./ipi_rna_seq_processor_base"

class Polyphemus::IpiRnaSeqMatrixProcessor < Polyphemus::IpiRnaSeqProcessorBase
  MATRIX_FILE_REGEX = /^gene_(counts|tpm)_table.tsv$/
  PROJECT = "ipi"
  BUCKET = "data"
  MAGMA_MODEL = "rna_seq"

  def initialize
    super
    @gene_ids_by_project = {}
  end

  def process(cursor, files)
    matrix_files = files.select do |file|
      file.file_name =~ MATRIX_FILE_REGEX
    end

    logger.info("Found #{matrix_files.length} matrix files: #{matrix_files.map { |f| f.file_path }.join(",")}")

    magma_gene_ids = matrix_gene_ids(cursor[:project_name])

    download_files(matrix_files) do |matrix_file, tmp_file|
      # The input matrices are transposed...
      # rows are gene_ids
      # columns are rna_seq tube_names.
      logger.info("Downloading #{matrix_file.file_path}.")
      csv = CSV.parse(::File.read(tmp_file.path), headers: true, col_sep: "\t")
      csv.by_col!
      attribute_name = matrix_file.file_name.sub("_table.tsv", "")

      data_gene_ids_map = Hash[csv[0].map.with_index.to_a]

      logger.info("Processing #{matrix_file.file_path}.")
      csv.each.with_index do |col, index|
        next if index == 0

        logger.info("Instantiating matrix class for index #{index.to_s}.")
        matrix = MagmaRnaSeqMatrix.new(
          raw_data: col,
          magma_gene_ids: magma_gene_ids,
          data_gene_ids_map: data_gene_ids_map,
          helper: @helper,
        )

        next if @helper.is_non_cancer_sample?(matrix.tube_name)

        logger.info("Preparing revision request for #{matrix.tube_name}.")
        update_for_cursor(cursor) do |update_request|
          update_request.update_revision(MAGMA_MODEL, matrix.tube_name, {
            "#{attribute_name}": matrix.to_array,
          })
          logger.info("Updating #{attribute_name} for record #{matrix.tube_name}.")
        end
        logger.info("Update complete.")
      end
    end

    logger.info("Done")
  end

  private

  def matrix_gene_ids(project_name)
    return @gene_ids_by_project[project_name] if @gene_ids_by_project.key?(project_name)

    @gene_ids_by_project[project_name] = magma_client.retrieve(Etna::Clients::Magma::RetrievalRequest.new(
      project_name: project_name,
      model_name: MAGMA_MODEL,
      attribute_names: [],
      record_names: [],
    )).models.model(MAGMA_MODEL).template.attributes.attribute("gene_counts").validation["value"]

    @gene_ids_by_project[project_name]
  end

  class MagmaRnaSeqMatrix
    def initialize(raw_data:, magma_gene_ids:, data_gene_ids_map:, helper:)
      @raw = raw_data
      @helper = helper
      @magma_gene_ids = magma_gene_ids
      @data_gene_ids_map = data_gene_ids_map
    end

    def raw_tube_name
      @raw[0]
    end

    def raw_gene_data
      @raw[1]
    end

    def tube_name
      return @helper.control_name(raw_tube_name) if @helper.is_control?(raw_tube_name)

      @helper.corrected_rna_seq_tube_name(raw_tube_name)
    end

    def to_array
      [].tap do |result|
        @magma_gene_ids.each do |magma_gene_id|
          value = @data_gene_ids_map.key?(magma_gene_id) ?
            raw_gene_data[@data_gene_ids_map[magma_gene_id]] :
            0

          result << value.to_f
        end
      end
    end
  end
end
