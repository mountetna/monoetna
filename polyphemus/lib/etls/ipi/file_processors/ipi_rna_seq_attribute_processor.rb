require_relative "./ipi_rna_seq_processor_base"

class Polyphemus::IpiRnaSeqAttributeProcessor < Polyphemus::IpiRnaSeqProcessorBase
  ATTRIBUTE_FILE_REGEX = /^rnaseq_table.tsv$/
  PROJECT = "ipi"
  BUCKET = "data"
  MAGMA_MODEL = "rna_seq"

  def process(cursor, files)
    attribute_files = files.select do |file|
      file.file_name =~ ATTRIBUTE_FILE_REGEX
    end

    logger.info("Found #{attribute_files.length} attribute files: #{attribute_files.map { |f| f.file_path }.join(",")}")

    return if attribute_files.empty?

    model_attributes = attributes(cursor[:project_name])

    download_files(attribute_files) do |attribute_file, tmp_file|
      tmp_file.readline # Get rid of the first, non-header row

      csv = CSV.new(tmp_file, headers: true, col_sep: "\t")

      csv.each do |row|
        rna_seq = MagmaRnaSeq.new(row, model_attributes)

        next if @helper.is_non_cancer_sample?(rna_seq.tube_name)

        update_for_cursor(cursor) do |update_request|
          update_request.update_revision(MAGMA_MODEL, rna_seq.tube_name, rna_seq.to_hash)
          logger.info("Updating record #{rna_seq.tube_name}.")
        end
      end
    end

    logger.info("Done")
  end

  private

  def attributes(project_name)
    magma_client.retrieve(Etna::Clients::Magma::RetrievalRequest.new(
      project_name: project_name,
      model_name: MAGMA_MODEL,
      attribute_names: [],
      record_names: [],
    )).models.model(MAGMA_MODEL).template.attributes
  end

  class MagmaRnaSeq
    # these values are too large for postgres ints
    ATTRIBUTES_TO_SKIP = [
      "raw_base_count", "filtered_base_count",
      "raw_read_count", "filtered_read_count",
      "filter_passing_bases", "aligned_bases",
    ]

    CONTROL_ATTRIBUTES_TO_SKIP = [
      "compartment",
    ]

    def initialize(table_row, attributes)
      @raw = table_row
      @helper = IpiHelper.new
      @attributes = attributes
    end

    def raw_tube_name
      @raw[0]
    end

    def tube_name
      return @helper.control_name(raw_tube_name) if @helper.is_control?(raw_tube_name)

      @helper.corrected_rna_seq_tube_name(raw_tube_name)
    end

    def method_missing(name, *args, &block)
      @raw[name.to_s]
    end

    def cell_number
      @raw["cell_count"]
    end

    def duplication_pct
      @raw["duplication_rate"]
    end

    def ribosomal_read_count
      @raw["reads"]
    end

    def input_read_count
      @raw["input_reads"]
    end

    def uniq_map_read_count
      @raw["uniq_map_reads"]
    end

    def multimap_lte20_read_count
      @raw["multimap_lte20_reads"]
    end

    def multimap_gt20_read_count
      @raw["multimapp_gt20_reads"]
    end

    def chimeric_read_count
      @raw["chimeric_reads"]
    end

    def chromosomal_read_count
      @raw["chromosomal"]
    end

    def mitochondrial_read_count
      @raw["mitochondrial"]
    end

    def median_5prime_bias
      @raw["median_5p_bias"]
    end

    def median_3prime_bias
      @raw["median_3p_bias"]
    end

    def eisenberg_score
      @raw["EHK"]
    end

    def expressed_eisenberg_genes
      @raw["expressed_EHK_genes"]
    end

    def raw_mean_length
      parts = @raw["raw_mean_length"].split(",")
      (parts.first.to_f + parts.last.to_f) / 2
    end

    def filtered_mean_length
      parts = @raw["filtered_mean_length"].split(",")
      (parts.first.to_f + parts.last.to_f) / 2
    end

    def compartment
      raw_value = @raw["compartment"]

      return raw_value if valid_compartment_values.include?(raw_value)

      # live2 => live
      return raw_value.match(compartment_increment_regex)[:compartment] if is_numeric_increment?(raw_value)

      # tpost, tpre, tpre2
      "other"
    end

    def valid_compartment_values
      @compartment_values ||= @attributes.attribute("compartment").validation["value"]
    end

    def compartment_increment_regex
      /^(?<compartment>(#{valid_compartment_values.join("|")}))\d+$/
    end

    def is_numeric_increment?(compartment)
      !!compartment.match(compartment_increment_regex)
    end

    def to_hash
      {}.tap do |result|
        @attributes.attribute_keys.each do |attr_name|
          attr = @attributes.attribute(attr_name)

          next if skip_attribute?(attr)

          value = self.send(attr_name.to_sym)

          value = value.to_i if attr.attribute_type == Etna::Clients::Magma::AttributeType::INTEGER
          value = value.to_f if attr.attribute_type == Etna::Clients::Magma::AttributeType::FLOAT

          result[attr_name.to_sym] = value
        end
      end.compact
    end

    def skip_attribute?(attr)
      attr.hidden ||
      attr.restricted ||
      attr.read_only ||
      attr.link_model_name ||
      ATTRIBUTES_TO_SKIP.include?(attr.attribute_name) ||
      (@helper.is_control?(raw_tube_name) && CONTROL_ATTRIBUTES_TO_SKIP.include?(attr.attribute_name))
    end
  end
end
