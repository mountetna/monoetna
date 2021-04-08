require_relative "magma_record_etl"
require_relative "hash_scan_based_etl_scanner"

class Polyphemus
  class MagmaRecordFileEtl < Polyphemus::MagmaRecordEtl
    def scanner
      HashScanBasedEtlScanner.new
    end

    def scanner_definition
      scanner_def = super

      scanner_def.start_batch_state do |cursor|
        # Use query instead of retrieve so we don't have to
        #   merge in a separate file md5 query later.
        query_request = Etna::Clients::Magma::QueryRequest.new(
          project_name: cursor[:project_name],
          query: [cursor[:model_name]],
          order: "updated_at",
          page: 1,
        )

        prepare_request(cursor, query_request)

        query_request
      end.execute_batch_find do |query_request, i|
        query_request.page_size = @limit * i

        begin
          self.magma_client.query(query_request).answer.map do |result|
            { result.first => result.last.map.with_index { |value, i|
              [queried_attribute_names[i], value]
            }.to_h }
          end
        rescue Etna::Error => e
          {} if e.message =~ /Page.*not found/
        end
      end.result_file_hashes do |record|
        record.values.first.slice(*file_attributes_to_query)
      end

      scanner_def
    end

    def prepare_request(cursor, request)
      fetch_template(cursor)

      # No time precision worries, so just use cursor's updated_at here
      request.query << ["updated_at", "::>=", (cursor.updated_at).iso8601] unless cursor.updated_at.nil?
      request.query << has_attribute_filters
      request.query << "::all"
      request.query << attribute_query_terms
    end

    def fetch_template(cursor)
      # Fetch the model template so that we can correctly
      #   identify any file attributes.
      retrieve_request = Etna::Clients::Magma::RetrievalRequest.new(
        project_name: cursor[:project_name],
        model_name: cursor[:model_name],
        attribute_names: "all",
        record_names: [],
      )
      @template = self.magma_client.retrieve(
        retrieve_request
      ).models.model(
        retrieve_request.model_name
      ).template
    end

    def has_attribute_filters
      ["::or"] + file_attributes_to_query.map { |a| ["::has", a] }
    end

    def attribute_query_terms
      return all_attribute_names if all_attributes?

      non_file_attribute_names + file_attributes_to_query.map { |a| [a, "::md5"] }
    end

    def queried_attribute_names
      @queried_attribute_names ||= attribute_query_terms.flatten.select { |a| a != "::md5" }
    end

    def non_file_attribute_names
      attribute_set = all_attributes? ? all_attribute_names : @attribute_names

      attribute_set - file_attributes_to_query
    end

    def all_attribute_names
      @all_attribute_names ||= @template.attributes.all.map { |a| a.name }
    end

    def all_file_attribute_names
      @all_file_attribute_names ||= all_attribute_names.select { |a|
        [Etna::Clients::Magma::AttributeType::FILE,
         Etna::Clients::Magma::AttributeType::IMAGE,
         Etna::Clients::Magma::AttributeType::FILE_COLLECTION].include?(a.attribute_type)
      }.map { |a| a.name }
    end

    def file_attributes_to_query
      @file_attributes_to_query ||= all_attributes? ?
        all_file_attribute_names :
        all_file_attribute_names.select { |a|
        @attribute_names.include?(a)
      }
    end

    def all_attributes?
      @attribute_names == "all"
    end
  end
end
