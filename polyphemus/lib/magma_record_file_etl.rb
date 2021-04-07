require_relative "magma_record_etl"
require_relative "time_scan_based_etl_scanner_with_hash"

class Polyphemus
  class MagmaRecordFileEtl < Polyphemus::MagmaRecordEtl
    def scanner
      TimeScanBasedEtlScannerWithHash.new
    end

    def scanner_definition
      scanner_def = super

      scanner_def.execute_batch_find do |retrieve_request, i|
        # Overwrite this to replace the File path / url
        #   data with the MD5 hash via Magma query

        retrieve_request.page_size = @limit * i
        magma_model = self.magma_client.retrieve(retrieve_request).models.model(retrieve_request.model_name)
        documents = magma_model.documents
        @template = magma_model.template

        raw_documents = documents.document_keys.map { |i| { i => documents.document(i) } }

        hashes = md5_hashes_for_records(
          retrieve_request.project_name,
          retrieve_request.model_name,
          raw_documents.map { |d| scanner_def.result_id(d) }
        )

        raw_documents.map do |doc|
          record_identifier = doc.keys.first
          doc[record_identifier].update(
            hashes.key?(record_identifier) ? hashes[record_identifier] : {}
          )
          doc
        end
      end.result_file_hashes do |record|
        record.slice(*(file_attribute_names.map { |n| n.to_sym }))
      end

      scanner_def
    end

    def md5_hashes_for_records(project_name, model_name, record_names)
      return {} unless file_attribute_names.length > 0

      query_request = Etna::Clients::Magma::QueryRequest.new(
        project_name: project_name,
        query: [model_name,
                ["::identifier", "::in", record_names],
                "::all",
                file_attribute_names.map { |a| [a, "::md5"] }],
      )

      self.magma_client.query(query_request).answer.map do |result|
        [result.first, result.last.map.with_index { |value, i|
          [file_attribute_names[i], value]
        }.to_h]
      end.to_h
    end

    def all_file_attribute_names
      @all_file_attribute_names ||= @template.attributes.all.select { |a|
        [Etna::Clients::Magma::AttributeType::FILE,
         Etna::Clients::Magma::AttributeType::IMAGE,
         Etna::Clients::Magma::AttributeType::FILE_COLLECTION].include?(a.attribute_type)
      }.map { |a| a.name }
    end

    def file_attribute_names
      @file_attribute_names ||= @attribute_names == "all" ?
        all_file_attribute_names :
        all_file_attribute_names.select { |a|
        @attribute_names.include?(a)
      }
    end
  end
end
