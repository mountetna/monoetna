# This command should probably be an ETL. We should convert if possible,
#   but I haven't done so because I can't tell what the "search and find
#   files" infrastructure should be.
# This command recursively searches a directory tree to find
#   rna_seq output files (that end with `_rna_seq_table.tsv`)
#   and puts those data into Magma records.

require 'csv'
require 'find'
require 'net/http'
require 'ostruct'
require_relative '../helpers'

class RnaSeqMetrics
  attr_reader :data, :magma_client, :rna_seq_model, :execute, :project_name
  def initialize(project_name, tsv_path, magma_client, rna_seq_model, execute=false)
    @tsv_path = tsv_path
    @project_name = project_name
    @magma_client = magma_client
    @rna_seq_model = rna_seq_model
    @execute = execute

    @data = {}
    process_tsv
  end

  def value(key)
    case key
    when :raw_mean_length
      data[key].split(',').first
    when :filtered_mean_length
      data[key].split(',').first
    when :rna_seq_plate
      "Plate#{data[key]}"
    when :cell_number
      data[:cell_count]
    when :median_5prime_bias
      data[:median_5p_bias]
    when :median_3prime_bias
      data[:median_3p_bias]
    when :eisenberg_score
      data[:EHK]
    when :expressed_eisenberg_genes
      data[:expressed_EHK_genes]
    when :uniq_map_read_count
      data[:uniq_map_reads]
    when :input_read_count
      data[:input_reads]
    when :multimap_lte20_read_count
      data[:multimap_lte20_reads]
    when :multimap_gt20_read_count
      data[:multimapp_gt20_reads]
    when :chimeric_read_count
      data[:chimeric_reads]
    when :chromosomal_read_count
      data[:chromosomal]
    when :mitochondrial_read_count
      data[:mitochondrial]
    when :flag
      data[:low_exp_flag] if data[:low_exp_flag].downcase != 'false'
    else
      data[key]
    end
  end

  def process_tsv
    CSV.foreach(@tsv_path, col_sep: "\t") do |tsv_line|
      data[:ribosomal_read_count] = tsv_line[2] if tsv_line[0] == 'ribosomal_rna' && tsv_line[1] == 'reads'
      data[tsv_line[1].to_sym] = tsv_line[2] unless tsv_line[1].nil?
    end
  end

  def filename
    @tsv_path.split('/').last
  end

  def is_control?
    filename =~ /control/i
  end

  def is_jurkat?
    is_control? && filename =~ /jurkat/i
  end

  def is_uhr?
    is_control? && filename =~ /uhr/i
  end

  def tube_name
    # Have to account for the validation, mostly for control?

    return "Control_Jurkat.Plate#{data[:rna_seq_plate]}" if is_jurkat?

    return "Control_UHR.Plate#{data[:rna_seq_plate]}" if is_uhr?

    "#{data[:sample].upcase}.rna.#{data[:compartment].downcase}"
  end

  def upload_rna_seq
    puts "Updating rna seq metrics: #{tube_name}"
    doc = create_update_doc
    puts doc
    if execute
      puts "Sending the update request."
      update_request = Etna::Clients::Magma::UpdateRequest.new(project_name: project_name)
      update_request.update_revision(
        'rna_seq',
        tube_name,
        doc)
      magma_client.update(update_request)
    end
  end

  def create_update_doc
    doc = {}
    uneditable_attributes = [
      'created_at', 'updated_at', 'tube_name',
      'gene_tpm', 'gene_counts', 'flag', 'unmapped_read_count',
      # These values are too large for Postgres integer field
      'aligned_bases', 'filter_passing_bases', 'filtered_base_count',
      'raw_base_count']
    uneditable_control_attributes = ['sample', 'compartment']
    rna_seq_model.template.attributes.attribute_keys.each do |attribute_name|
      next if uneditable_attributes.include?(attribute_name)
      next if is_control? && uneditable_control_attributes.include?(attribute_name)

      doc[attribute_name] = value(attribute_name.to_sym)
    end

    doc
  end
end

class RnaSeqGenes
  attr_reader :data, :magma_client, :rna_seq_model, :execute, :project_name
  def initialize(project_name, tsv_path, magma_client, rna_seq_model, execute=false)
    @tsv_path = tsv_path
    @project_name = project_name
    @magma_client = magma_client
    @rna_seq_model = rna_seq_model
    @execute = execute

    @data = {}
    process_tsv
  end

  def filename
    @tsv_path.split('/').last
  end

  def is_control?
    filename =~ /control/i
  end

  def is_jurkat?
    is_control? && filename =~ /jurkat/i
  end

  def is_uhr?
    is_control? && filename =~ /uhr/i
  end

  def tube_name
    # Have to account for the validation, mostly for control?

    from_filename = filename.sub('.rsem.genes.results', '')

    return from_filename unless is_jurkat? || is_uhr?

    plate_number = /.*plate(?<plate_number>\d+)/i.match(from_filename).named_captures['plate_number']

    return "Control_Jurkat.Plate#{plate_number}" if is_jurkat?

    return "Control_UHR.Plate#{plate_number}" if is_uhr?

    from_filename
  end

  def matrix_columns
    # assume same columns for each gene attribute
    @matrix_columns ||= rna_seq_model.template.attributes.attribute('gene_tpm').options
  end

  def process_tsv
    CSV.foreach(@tsv_path, col_sep: "\t", headers: true) do |tsv_line|
      data[tsv_line["gene_id"]] = {
        gene_counts: tsv_line["expected_count"],
        gene_tpm: tsv_line["TPM"]
      }
    end
  end

  def get_data_as_array(attribute)
    # The IPI gene_tpm and gene_counts validation options have
    #   an outdated Ensembl gene set, which includes 280551
    #   genes. However our dataset only includes 27992 genes.
    # Here we insert `0` for those extra 59 genes.
    # NOTE: We can't take them out of the validation without
    #   risk of ruining the existing data.
    array = []
    matrix_columns.each do |gene_id|
      if data.key?(gene_id)
        array << data[gene_id][attribute].to_f
      else
        array << 0
      end
    end

    array
  end

  def tpms
    get_data_as_array(:gene_tpm)
  end

  def counts
    get_data_as_array(:gene_counts)
  end

  def upload_rna_seq
    puts "Updating rna seq genes: #{tube_name}"

    doc = {
      gene_counts: counts,
      gene_tpm: tpms
    }

    # puts doc
    if execute
      puts "Sending the update request."
      update_request = Etna::Clients::Magma::UpdateRequest.new(project_name: project_name)
      update_request.update_revision(
        'rna_seq',
        tube_name,
        doc)
      magma_client.update_json(update_request)
    end
  end
end

class ProcessRnaSeqOutput < Struct.new(:magma_client, :project_name, :file_path, :execute, keyword_init: true)
  include WithLogger

  def initialize(**opts)
    super(**{project_name: 'ipi', execute: false}.update(opts))
  end

  def rna_seq_model
    @rna_seq_model ||= begin
      params = {
        model_name: 'rna_seq',
        record_names: [],
        attribute_names: 'all'
      }

      request = Etna::Clients::Magma::RetrievalRequest.new(
        project_name: project_name, **params)
      magma_client.retrieve(request).models.model('rna_seq')
    end
  end

  def load_rna_seq_metrics(filepath)
    puts "Processing metrics file #{filepath}."
    rna_seq_metrics = RnaSeqMetrics.new(
      project_name, filepath, magma_client, rna_seq_model, execute)
    rna_seq_metrics.upload_rna_seq()
  end

  def load_rna_seq_genes(filepath)
    puts "Processing genes file #{filepath}."
    rna_seq_genes = RnaSeqGenes.new(
      project_name, filepath, magma_client, rna_seq_model, execute)
    rna_seq_genes.upload_rna_seq()
  end

  def process_rna_seq_metrics
    Find.find(file_path) do |path|
      load_rna_seq_metrics(path) if path =~ /.*_rnaseq_table\.tsv$/
    end
  end

  def process_rna_seq_genes
    Find.find(file_path) do |path|
      load_rna_seq_genes(path) if path =~ /.*\.rsem\.genes\.results$/
    end
  end

  def process_rna_seq
    process_rna_seq_metrics
    process_rna_seq_genes
  end
end