require_relative "../project_watch_folders_etl"
require_relative "../project_watch_files_etl"
require_relative "./linkers/mvir1_rna_seq_raw_fastq_files_linker"
require_relative "./record_creators/mvir1_rna_seq_record_creator"
require_relative "../linker_processor"

class Polyphemus
  module Mvir
    class WatchFoldersConfig < Polyphemus::ProjectWatchFoldersConfig
      def initialize
        super(project_name: "mvir1")

        process_folders_with(
          bucket("data").watcher("create_rna_seq_records").watch(/^bulk_RNASeq\/raw\/.*$/),
          Polyphemus::MvirRnaSeqRecordCreator.new
        )

        process_watch_type_with(
          bucket("data").watcher("link_rna_seq_raw_fastq_files").watch(/^bulk_RNASeq\/raw\/.*$/),
          Polyphemus::LinkerProcessor.new(linker: rna_seq_fastq_linker, model_name: "rna_seq")
        )
      end

      def rna_seq_fastq_linker
        Polyphemus::MvirRnaSeqRawFastqFilesLinker.new(project_name: project_name, bucket_name: "data")
      end
    end

    class MvirWatchFoldersEtl < Polyphemus::ProjectWatchFoldersEtl
      def initialize
        super(WatchFoldersConfig.new, limit: 40)
      end
    end

    class MvirWatchFilesEtl < Polyphemus::ProjectWatchFilesEtl
      def initialize
        super(WatchFoldersConfig.new)
      end
    end

    class MvirPropagateFolderUpdatedAt < Polyphemus::ProjectPropagateFolderUpdatedAtEtl
      def initialize
        super(WatchFoldersConfig.new, limit: 40)
      end
    end
  end
end
