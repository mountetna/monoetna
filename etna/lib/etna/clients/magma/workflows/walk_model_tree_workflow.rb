require 'ostruct'
require 'digest'
require 'fileutils'
require 'tempfile'

module Etna
  module Clients
    class Magma
      class WalkModelTreeWorkflow < Struct.new(:magma_crud, :logger, keyword_init: true)
        def initialize(**args)
          super(**({}.update(args)))
          @template_for = {}
        end

        def masked_attributes(template:, model_attributes_mask:, model_name:)
          attributes_mask = model_attributes_mask[model_name]
          return ["all", "all"] if attributes_mask.nil?
          [(attributes_mask + [template.identifier, 'parent']).uniq, attributes_mask]
        end

        def attribute_included?(mask, attribute_name)
          return true if mask == "all"
          mask.include?(attribute_name)
        end

        def template_for(model_name)
          @template_for[model_name] ||= magma_crud.magma_client.retrieve(RetrievalRequest.new(
              project_name: magma_crud.project_name,
              model_name: model_name,
              record_names: [],
              attribute_names: [],
          )).models.model(model_name).template
        end

        def walk_from(
            model_name,
            record_names = 'all',
            model_attributes_mask: {},
            model_filters: {},
            page_size: 100,
            &block)
          q = [ { model_name: model_name, from: nil, record_names: record_names } ]
          seen = Set.new

          while (path = q.pop)
            model_name = path[:model_name]
            next if seen.include?([path[:from], model_name])
            seen.add([path[:from], model_name])

            template = template_for(model_name)
            query_attributes, walk_attributes = masked_attributes(template: template, model_attributes_mask: model_attributes_mask, model_name: model_name)

            request = RetrievalRequest.new(
                project_name: magma_crud.project_name,
                model_name: model_name,
                record_names: path[:record_names],
                filter: model_filters[model_name],
                attribute_names: query_attributes,
                page_size: page_size, page: 1
            )

            related_models = {}

            magma_crud.page_records(model_name, request) do |response|
              tables = []
              collections = []
              links = []
              attributes = []

              model = response.models.model(model_name)

              template.attributes.attribute_keys.each do |attr_name|
                next unless attribute_included?(query_attributes, attr_name)
                attributes << attr_name

                attr = template.attributes.attribute(attr_name)
                if attr.attribute_type == AttributeType::TABLE
                  tables << attr_name
                elsif attr.attribute_type == AttributeType::COLLECTION
                  related_models[attr.link_model_name] ||= Set.new
                  collections << attr_name
                elsif attr.attribute_type == AttributeType::LINK
                  related_models[attr.link_model_name] ||= Set.new
                  links << attr_name
                elsif attr.attribute_type == AttributeType::CHILD
                  related_models[attr.link_model_name] ||= Set.new
                  links << attr_name
                elsif attr.attribute_type == AttributeType::PARENT && attribute_included?(walk_attributes, attr_name)
                  related_models[attr.link_model_name] ||= Set.new
                  links << attr_name
                end
              end

              model.documents.document_keys.each do |key|
                record = model.documents.document(key).slice(*attributes)

                # Inline tables inside the record
                tables.each do |table_attr|
                  record[table_attr] = record[table_attr].map do |id|
                    response.models.model(template.attributes.attribute(table_attr).link_model_name).documents.document(id)
                  end unless record[table_attr].nil?
                end

                collections.each do |collection_attr|
                  record[collection_attr].each do |collected_id|
                    related_models[template.attributes.attribute(collection_attr).link_model_name].add(collected_id)
                  end unless record[collection_attr].nil?
                end

                links.each do |link_attr|
                  related_models[template.attributes.attribute(link_attr).link_model_name].add(record[link_attr]) unless record[link_attr].nil?
                end

                yield template, record
              end
            end

            related_models.each do |link_model_name, id_set|
              next if id_set.empty?
              q.push({ model_name: link_model_name, from: model_name, record_names: id_set.to_a })
            end
          end
        end
      end
    end
  end
end

