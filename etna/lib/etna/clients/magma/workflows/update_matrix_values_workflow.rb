module Etna
  module Clients
    class Magma
      class UpdateMatrixValuesWorkflow < Struct.new(:magma_client, :project_name, :filepath, :model_name, :attribute_name, :execute, :logger, :tsv, keyword_init: true)
        def initialize(opts)
          super(**{}.update(opts))
        end

        def model
          @model ||= begin
              params = {
                model_name: self.model_name,
                record_names: [],
                attribute_names: "all",
              }

              request = Etna::Clients::Magma::RetrievalRequest.new(
                project_name: project_name, **params,
              )
              self.magma_client.retrieve(request).models.model(self.model_name)
            end
        end

        def filename
          ::File.basename(self.filepath)
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
          # Could also take this from the header for column 2...

          from_filename = filename.sub("_counts.#{tsv ? "tsv" : "csv"}", "")

          return from_filename unless is_jurkat? || is_uhr?

          plate_number = /.*plate(?<plate_number>\d+)/i.match(from_filename).named_captures["plate_number"]

          return "Control_Jurkat.Plate#{plate_number}" if is_jurkat?

          return "Control_UHR.Plate#{plate_number}" if is_uhr?

          from_filename
        end

        def matrix_columns
          # assume same columns for each gene attribute
          @matrix_columns ||= model.template.attributes.attribute(self.attribute_name).validation["value"]
        end

        def data
          @data ||= begin
              @data = {}
              CSV.foreach(self.filepath, col_sep: tsv ? "\t" : ",", headers: true) do |csv_line|
                @data[csv_line[0]] = csv_line[1]
              end

              @data
            end
        end

        def get_data_as_array
          # The gene_tpm and gene_counts validation options have
          #   an outdated Ensembl gene set, which includes 280551
          #   genes. However our dataset only includes 27992 genes.
          # Here we insert `0` for those extra 59 genes.
          # NOTE: We can't take them out of the validation without
          #   risk of ruining the existing data.
          array = []
          matrix_columns.each do |gene_id|
            if data.key?(gene_id)
              array << data[gene_id].to_f
            else
              array << 0
            end
          end

          array
        end

        def upload_values
          puts "Updating #{self.model_name} #{self.attribute_name}: #{tube_name}"

          doc = {
            "#{self.attribute_name}": get_data_as_array,
          }

          logger&.debug(doc)

          if execute
            puts "Sending the update request."
            update_request = Etna::Clients::Magma::UpdateRequest.new(project_name: self.project_name)
            update_request.update_revision(
              self.model_name,
              tube_name,
              doc
            )
            self.magma_client.update_json(update_request)
          end
        end
      end
    end
  end
end
