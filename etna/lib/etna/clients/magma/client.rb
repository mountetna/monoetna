require 'net/http/post/multipart'
require 'singleton'
require_relative '../base_client'
require_relative '../../client'
require_relative './models'

module Etna
  module Clients
    class LocalMagmaClient # May only be used from within magma app.
      def logger
        ::Magma.instance.logger
      end

      def user
        permissions = Etna::Permissions.new([])
        permissions.add_permission(Etna::Permission.new('A', 'administration'))
        @user ||= Etna::User.new({name: 'Admin', email: 'etnaagent@ucsf.edu', perm: permissions.to_string})
      end

      def fake_request(request)
        Rack::Request.new({
          'rack.request.params' => request.as_json,
          'etna.logger' => self.logger,
          'etna.user' => user,
          'etna.request_id' => 'local',
        })
      end

      def parse_json_response(response)
        status, headers, body = response
        body = body.join('')
        if status < 200 || status >= 300
          raise Etna::Error.new(body, status)
        end

        JSON.parse(body)
      end

      def retrieve(retrieval_request = RetrievalRequest.new)
        controller = ::RetrieveController.new(fake_request(retrieval_request), nil)
        Magma::RetrievalResponse.new(parse_json_response(controller.action))
      end

      def update_model(update_model_request = UpdateModelRequest.new)
        controller = ::UpdateModelController.new(fake_request(update_model_request), nil)
        Magma::UpdateModelResponse.new(parse_json_response(controller.action))
      end
    end

    class Magma < Etna::Clients::BaseClient
      # This endpoint returns models and records by name:
      # e.g. params:
      # {
      #   model_name: "model_one", # or "all"
      #   record_names: [ "rn1", "rn2" ], # or "all",
      #   attribute_names:  "all"
      # }
      def retrieve(retrieval_request = RetrievalRequest.new)
        json = nil
        @etna_client.post('/retrieve', retrieval_request) do |res|
          json = JSON.parse(res.body)
        end

        RetrievalResponse.new(json)
      end

      # This 'query' end point is used to fetch data by graph query
      # See question.rb for more detail
      def query(query_request = QueryRequest.new)
        json = nil
        @etna_client.post('/query', query_request) do |res|
          json = JSON.parse(res.body)
        end

        QueryResponse.new(json)
      end

      def update(update_request = UpdateRequest.new)
        json = nil

        @etna_client.multipart_post('/update', update_request.encode_multipart_content) do |res|
          json = JSON.parse(res.body)
        end

        UpdateResponse.new(json)
      end

      class Pager
        def initialize(client, update_request, page_size)
          @client = client

          @update_request = update_request

          @page_size = [ page_size || 2, 2 ].max
        end

        def response
          page_count = 0
          page = {}

          results = UpdateResponse.new({"models" => {}})

          flat_records.each do |flat_record|
            page_count += 1 + flat_record[:tables].length

            add_to_page( page, flat_record )

            flat_record[:tables].each do |table_record|
              add_to_page(page, table_record)
            end

            if page_count > @page_size
              results.update(update_part(page))
              page = {}
              page_count = 0
            end
          end

          results.update(update_part(page))

          return results
        end

        def revisions
          @update_request.revisions
        end

        def flat_records
          project = @client.retrieve({
            project_name: @update_request.project_name,
            hide_templates: false
          })

          flat_records = {}
          revisions.each do |model_name, records|
            records.each do |record_name, record|
              flat = {
                  model_name: model_name,
                  record_name: record_name,
                  record: record
              }
              flat_name = "#{model_name}.#{record_name}"

              # if the update sets a table attribute, attach table records to their parent so they are submitted together
              flat_tables = []
              record.each do |attribute_name, value|
                att = project.models.model(model_name).template.attributes.attribute(attribute_name)
                if att.attribute_type == 'table'
                  value.each do |table_entry_name|
                    table_record = revisions[ att.link_model_name ][table_entry_name]
                    flat_tables.push(
                      {
                        model_name: att.link_model_name,
                        record_name: table_entry_name,
                        record: table_record
                      }
                    )

                    # exclude the attached table record from the flattened list of records
                    flat_records[ "#{ att.link_model_name }.#{table_entry_name}" ] = nil
                  end
                end
              end
              flat[:tables] = flat_tables
              flat_records[ flat_name ] = flat unless flat_records.has_key?(flat_name)
            end
          end

          return flat_records.values.compact.shuffle
        end

        def update_part(page)
          @client.update_json(
            UpdateRequest.new(
              project_name: @update_request.project_name,
              revisions: page,
              dry_run: @update_request.dry_run,
              autolink: @update_request.autolink
            )
          )
        end

        def add_to_page(page, entry)
          page[ entry[:model_name] ] ||= {}
          page[ entry[:model_name] ][ entry[:record_name] ] = entry[:record]
        end
      end

      def update_json(update_request = UpdateRequest.new, page_size=nil)
        json = nil

        if page_size
          return Etna::Clients::Magma::Pager.new(self, update_request, page_size).response
        end

        @etna_client.post('/update', update_request) do |res|
          json = JSON.parse(res.body)
        end

        UpdateResponse.new(json)
      end

      def update_model(update_model_request = UpdateModelRequest.new)
        json = nil
        @etna_client.post('/update_model', update_model_request) do |res|
          json = JSON.parse(res.body)
        end

        UpdateModelResponse.new(json)
      end
    end
  end
end
