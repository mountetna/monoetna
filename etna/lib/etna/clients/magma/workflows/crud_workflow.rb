require 'ostruct'

module Etna
  module Clients
    class Magma
      class MagmaCrudWorkflow < Struct.new(:magma_client, :project_name, :read_only, :logger, :batch_size, :cache_policy, keyword_init: true)
        attr_reader :recorded_updates

        def initialize(args)
          super(**{logger: Etna::Logger.new('/dev/stdout', 999999, 1024 * 1024), batch_size: 30, cache_policy: :update_on_write}.update(args))
        end

        def lookup_record(model_name, record_id)
          if (cached = (@cache ||= {}).dig(model_name, record_id))
            return cached
          end

          result = magma_client.retrieve(RetrievalRequest.new(project_name: project_name, record_names: [record_id], model_name: model_name))\
            .models.model(model_name).documents.document(record_id)

          if cache_policy
            ((@cache ||= {})[model_name] ||= {})[record_id = result]
          end

          result
        end

        def page_records(model_name, request = Etna::Clients::Magma::RetrievalRequest.new(
                project_name: project_name,
                model_name: model_name,
                record_names: 'all',
                page_size: 20,
                page: 1,
            ), &block)

          documents = Documents.new({})
          last_page = nil
          while last_page.nil? || last_page.models.model_keys.map { |k| last_page.models.model(k).documents.raw.length }.sum > 0
            attempts = 0
            begin
              attempts += 1
              last_page = magma_client.retrieve(request)
            # Unfortunately, paging in magma is not great and times out from time to time.
            rescue Net::ReadTimeout => e
              if attempts > 5
                raise e
              end

              retry
            rescue Etna::Error => e
              if e.status === 502
                if attempts > 5
                  raise e
                end

                retry
              else
                raise e unless e.message.include?('not found')
                break
              end
            end

            documents += last_page.models.model(model_name).documents unless block_given?
            yield last_page if block_given?
            request.page += 1
          end

          documents
        end

        # Todo: Introduce associative concatenation operations for response objects and return
        # one response that munges the batched responses together.
        def update_records(method: :update)
          @recorded_updates ||= UpdateRequest.new(project_name: project_name)

          request = UpdateRequest.new(project_name: project_name)
          yield request

          case cache_policy
          when :update_on_write
            @cache ||= {}
            request.revisions.each do |model_name, revisions|
              model_revisions = @cache[model_name] ||= {}

              revisions.each do |record_name, revision|
                if model_revisions.include? record_name
                  model_revisions[record_name].update(revision)
                end
              end
            end
          end

          revisions = request.revisions

          responses = []
          revisions.to_a.each_slice(batch_size) do |batch|
            request.revisions = batch.to_h
            magma_client.send(method, request) unless read_only
            responses << @recorded_updates.revisions.update(request.revisions)
          end

          responses
        end
      end
    end
  end
end

