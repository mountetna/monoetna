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

            request = RetrievalRequest.new(
                project_name: magma_crud.project_name,
                model_name: model_name,
                record_names: path[:record_names],
                filter: model_filters[model_name],
                page_size: page_size, page: 1
            )

            related_models = {}

            magma_crud.page_records(model_name, request) do |response|
              model = response.models.model(model_name)
              template = model.template

              tables = []
              collections = []
              links = []
              attributes = []

              template.attributes.attribute_keys.each do |attr_name|
                attributes_mask = model_attributes_mask[model_name]
                black_listed = !attributes_mask.nil? && !attributes_mask.include?(attr_name)
                next if black_listed && attr_name != template.identifier && attr_name != 'parent'
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
                elsif attr.attribute_type == AttributeType::PARENT && !black_listed
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

