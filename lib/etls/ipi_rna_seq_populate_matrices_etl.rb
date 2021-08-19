require "csv"
require_relative "../metis_file_for_magma_model_etl"
require_relative "../ipi/ipi_helper"

class Polyphemus::IpiRnaSeqPopulateMatricesEtl < Polyphemus::MetisFileForMagmaModelEtl
  # For getting the results back, we'll call the "plate number" the record_name
  PATH_REGEX = /.*\/(?<record_name>.*)_rnaseq_new\/results\/gene_(counts|tpm)_table.tsv$/
  PROJECT = "ipi"
  BUCKET = "data"
  CURSOR_MODEL = "rna_seq_plate"
  MAGMA_MODEL = "rna_seq"

  def initialize
    @helper = IpiHelper.new
    @gene_ids_by_project = {}

    super(
      project_bucket_model_tuples: [[PROJECT, BUCKET, CURSOR_MODEL]],
      file_name_globs: ["bulkRNASeq/**/gene_*_table.tsv"],
      metis_path_to_record_name_regex: PATH_REGEX,
      record_name_gsub_pair: [/plate/, "Plate"],
    )
  end

  def process(cursor, result_files_by_plate)
    logger.info("Processing plates: #{result_files_by_plate.map { |f| f.record_name }.join(",")}")

    update_request = Etna::Clients::Magma::UpdateRequest.new(
      project_name: cursor[:project_name],
    )

    result_files_by_plate.each do |plate_wrapper|
      plate_wrapper.files.each do |matrix_file|
        Tempfile.create do |tmp|
          metis_client.download_file(matrix_file) do |chunk|
            tmp << chunk
          end

          tmp.rewind

          # The input matrices are transposed...
          # rows are gene_ids
          # columns are rna_seq tube_names.
          csv = CSV.parse(::File.read(tmp.path), headers: true, col_sep: "\t")
          csv.by_col!
          attribute_name = matrix_file.file_name.sub("_table.tsv", "")

          data_gene_ids = csv[0]

          csv.each.with_index do |col, index|
            next if index == 0

            matrix = MagmaRnaSeqMatrix.new(
              raw_data: col,
              plate_name: plate_wrapper.record_name,
              magma_gene_ids: matrix_gene_ids(cursor[:project_name]),
              data_gene_ids: data_gene_ids,
            )

            next if @helper.is_non_cancer_sample?(matrix.tube_name)

            update_request.update_revision(MAGMA_MODEL, matrix.tube_name, {
              "#{attribute_name}": matrix.to_array,
            })
            logger.info("Updating #{attribute_name} for record #{matrix.tube_name}.")
          end
        end
      end
    end

    magma_client.update_json(update_request)

    logger.info("Done")
  end

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
    def initialize(raw_data:, plate_name:, magma_gene_ids:, data_gene_ids:)
      @raw = raw_data
      @helper = IpiHelper.new
      @magma_gene_ids = magma_gene_ids
      @plate_name = plate_name
      @data_gene_ids = data_gene_ids
    end

    def raw_tube_name
      @raw[0]
    end

    def raw_gene_data
      @raw[1]
    end

    def tube_name
      return @helper.control_name(raw_tube_name) if @helper.is_control?(raw_tube_name)

      @helper.corrected_rna_seq_tube_name(@plate_name, raw_tube_name)
    end

    def to_array
      [].tap do |result|
        @magma_gene_ids.each do |magma_gene_id|
          data_index = @data_gene_ids.find_index(magma_gene_id)

          value = data_index.nil? ? 0 : raw_gene_data[data_index]

          result << value.to_f
        end
      end
    end
  end
end
