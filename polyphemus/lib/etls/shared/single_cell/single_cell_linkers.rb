class Polyphemus
  class SingleCellRawLinker < Polyphemus::MetisFilesLinkerBase
  end

  class SingleCellProcessedLinker < Polyphemus::MetisFilesLinkerBase
    # Record name is two directories above a file
    RECORD_NAME_REGEX = /.*\/(?<record_name>.*)\/[^\/]*\/[^\/]*$/

    attr_reader :record_name_regex

    def initialize(
      project_name:,
      bucket_name:,
      attribute_regex_overrides: {},
      record_name_regex: RECORD_NAME_REGEX
    )
      super(project_name: project_name, bucket_name: bucket_name)
      @attribute_regex_overrides = attribute_regex_overrides
      @record_name_regex = record_name_regex
    end

    def link(model_name:, files:, cursor: nil)
      super(
        model_name: model_name,
        files_by_record_name: organize_metis_files_by_magma_record(
          metis_files: files,
          magma_record_names: current_magma_record_names(project_name, model_name),
          path_regex: @record_name_regex,
        ),
        attribute_regex: attribute_regex,
        cursor: cursor,
      )
    end

    def attribute_regex
      {
        # # only if CITE-seq
        "cite_antibody_key": /^ADT_keep_features\.list$/,
        "mux_index_key": /^IDX_map\.tsv$/,
        "processed_robject_rdata": /.*_scTransformed_processed\.RData$/,
        "processed_umap": /.*_umap\.pdf$/,
        "processing_pipeline_parameters": /^cutoffs\.yml$/,
        "filtered_counts_h5": /^filtered_feature_bc_matrix\.h5$/,
        "raw_counts_h5": /^raw_feature_bc_matrix\.h5$/,
        "tenx_aligned_bam": /^possorted_genome_bam\.bam$/,
        "tenx_aligned_bam_index": /^possorted_genome_bam\.bam\.bai$/,
        "tenx_cloupe_file": /^cloupe\.cloupe$/,
        "tenx_metrics_csv": /^metrics_summary\.csv$/,
        "tenx_molecule_info_h5": /^molecule_info\.h5$/,
        "tenx_web_summary": /^web_summary\.html$/,
      }.update(@attribute_regex_overrides)
    end

    def corrected_record_name(record_name)
    end

    def should_skip_record?(record_name)
      false
    end
  end
end