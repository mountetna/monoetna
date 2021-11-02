require_relative '../project_watch_folders_etl'
require_relative '../project_watch_files_etl'
require_relative '../project_propagate_folder_updated_at_etl'
require_relative './file_processors/ipi_rna_seq_attribute_processor'
require_relative './file_processors/ipi_rna_seq_matrix_processor'
require_relative './linkers/ipi_rna_seq_processed_files_linker'
require_relative './linkers/ipi_rna_seq_raw_fastq_files_linker'
require_relative './record_creators/ipi_rna_seq_and_plate_record_creator'
require_relative '../linker_processor'
require_relative '../shared/single_cell/single_cell_linkers'

class Polyphemus
  module Ipi
    module SingleCellLinkers
      def add_single_cell_linkers!
        process_folders_with(
          process_watch_type_with(
            bucket('data')
              .watcher('single_cell_pool_processed')
              .watch(/^single_cell_[^\/]*\/processed\/.*\/[^\/]+POOL[^\/]+\/.*$/),
            Polyphemus::LinkerProcessor.new(
              linker: single_cell_pooled_linker,
              model_name: 'sc_rna_seq_pool'
            ),
          ),
          SingleCellFolderProcessor.new(
            model: 'sc_rna_seq_pool',
            linker: single_cell_pooled_linker
          )
        )

        process_folders_with(
          process_watch_type_with(
            bucket('data')
              .watcher('single_cell_processed')
              .watch(/^single_cell_[^\/]*\/processed\/((?!POOL).)*$/),
            Polyphemus::LinkerProcessor.new(
              linker: single_cell_non_pooled_linker,
              model_name: 'sc_rna_seq',
            ),
          ),
          SingleCellFolderProcessor.new(
            model: 'sc_rna_seq',
            linker: single_cell_non_pooled_linker
          )
        )
      end

      def single_cell_non_pooled_linker
        SingleCellLinker.new(
          record_name_regex: /.*\/(?<record_name>IPI[^\/]*)\/.*$/
        )
      end

      def single_cell_pooled_linker
        SingleCellLinker.new
      end

      class SingleCellFolderProcessor
        include WithEtnaClients
        include WithSlackNotifications
        include WithLogger

        def initialize(
          model:,
          linker: SingleCellLinker.new
        )
          @model = model
          @linker = linker
        end

        VERSION_DIGIT_REGEX = /.*\/(v2_chemistry|v3_chemistry)\/.*/
        PRIME_DIGIT_REGEX = /.*\/(3prime_|5prime_).*/

        def calculate_record(folder, record_name)
          biospecimen = record_name.split('.').last
          prime_digit = nil
          version_digit = '3'

          match = VERSION_DIGIT_REGEX.match(folder.folder_path)
          if match
            if match[1].include?('2')
              version_digit = '2'
            elsif match[1].include?('1')
              version_digit = '1'
            else
              version_digit = '3'
            end
          end

          match = PRIME_DIGIT_REGEX.match(folder.folder_path)
          if match
            if match[1].include?('3')
              prime_digit = '3'
            else
              prime_digit = '5'
            end
          end

          if prime_digit.nil?
            logger.warn("Skipping #{folder.folder_path}, could not find prime component")
            return nil
          end

          chemistry = "10X_#{prime_digit}prime_v#{version_digit}"

          {
            'biospecimen' => biospecimen,
            'chemistry' => chemistry
          }
        end

        def process(cursor, folders)
          record_names = @linker.current_magma_record_names(cursor[:project_name], @model)

          matched_folders = folders.map do |folder|
            record_name = @linker.record_name_by_path(
              folder.folder_path + '/',
              @linker.record_name_regex
            )
            [folder, record_name] if record_names.include?(record_name)
          end.select { |f| f }

          matched_folders.each do |folder, record_name|
            record = calculate_record(folder, record_name)
            next if record.nil?

            begin
              magma_client.update_json(Etna::Clients::Magma::UpdateRequest.new(revisions: {
                @model => {
                  record_name => record
                }
              }, project_name: cursor[:project_name]))
            rescue Exception => e
              notify_slack(
                "Error processing ipi record directory #{folder.folder_path}.\n#{e.message}.",
                channel: "data-ingest-errors",
              )
              logger.log_error(e)
            end
          end
        end
      end

      class SingleCellLinker < Polyphemus::SingleCellProcessedLinker
        def initialize(**kwds)
          super(
            project_name: 'ipi',
            bucket_name: 'data',
            **kwds
          )
        end

        def corrected_record_name(record_name)
          if record_name.include?('POOL')
            record_name.gsub(/_/, '.').sub(/^([^.]*\.[^.]*\.)(.*)$/) do
              first = $1
              second = $2
              first + second.downcase
            end
          else
            record_name.sub(/^([^_.]*[_.][^_.]*[_.])([^_.]*)(.*)$/) do
              first = $1
              second = $2
              third = $3
              first.gsub(/_/, '.') + second.downcase + third
            end
          end
        end
      end
    end

    class WatchFoldersConfig < Polyphemus::ProjectWatchFoldersConfig
      include SingleCellLinkers

      def initialize
        super(project_name: 'ipi')

        add_single_cell_linkers!

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
          bucket('integral_data').watcher('link_rna_seq_raw_fastq_files').watch(/^(.*\/)?BulkRNASeq\/.*$/),
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