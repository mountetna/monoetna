require_relative '../project_watch_folders_etl'
require_relative '../project_watch_files_etl'
require_relative '../project_propagate_folder_updated_at_etl'
require_relative './file_processors/ipi_rna_seq_attribute_processor'
require_relative './file_processors/ipi_rna_seq_matrix_processor'
require_relative './linkers/ipi_rna_seq_processed_files_linker'
require_relative './linkers/ipi_rna_seq_raw_fastq_files_linker'
require_relative './record_creators/ipi_rna_seq_and_plate_record_creator'
require_relative '../linker_processor'

class Polyphemus
  module Ipi
    class WatchFoldersConfig < Polyphemus::ProjectWatchFoldersConfig
      def initialize
        super(project_name: 'ipi')

        process_watch_type_with(
          bucket('data').watcher('process_bulk_rna_seq_results').watch(/^bulkRNASeq\/.*\/results$/),
          Polyphemus::IpiRnaSeqAttributeProcessor.new,
            Polyphemus::IpiRnaSeqMatrixProcessor.new,
        )

        process_folders_with(
          process_watch_type_with(
            bucket('data').watcher('link_processed_rna_seq_files').watch(/^bulkRNASeq\/.*\/output\/.*$/),
            Polyphemus::LinkerProcessor.new(linker: rna_seq_processed_linker, model_name: 'rna_seq'),
          ),
          Polyphemus::IpiRnaSeqAndPlateRecordCreator.new
        )

        process_watch_type_with(
          bucket('integral_data').watcher('link_rna_seq_raw_fastq_files').watch(/^.*\/BulkRNASeq\/.*$/),
          Polyphemus::LinkerProcessor.new(linker: rna_seq_fastq_linker, model_name: 'rna_seq')
        )
      end

      def rna_seq_processed_linker
        Polyphemus::IpiRnaSeqProcessedFilesLinker.new(project_name: project_name, bucket_name: 'data')
      end

      def rna_seq_fastq_linker
        Polyphemus::IpiRnaSeqRawFastqFilesLinker.new(project_name: project_name, bucket_name: 'integral_data')
      end

      def ipi_rna_seq_and_plate_record_creator
        Polyphemus::IpiRnaSeqAndPlateRecordCreator.new
      end
    end

    class IpiWatchFoldersEtl < Polyphemus::ProjectWatchFoldersEtl
      def initialize
        super(WatchFoldersConfig.new)
      end
    end

    class IpiWatchFilesEtl < Polyphemus::ProjectWatchFilesEtl
      def initialize
        super(WatchFoldersConfig.new)
      end
    end

    class IpiPropagateFolderUpdatedAt < Polyphemus::ProjectPropagateFolderUpdatedAtEtl
      def initialize
        super(WatchFoldersConfig.new)
      end
    end
  end
end