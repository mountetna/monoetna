require 'ostruct'

module Etna
  module Clients
    class Magma
      class MagmaCrudWorkflow < Struct.new(:magma_client)
        def lookup_record(model_name, record_id)
          magma_client.retrieve(RetrievalRequest.new(project_name: project_name, record_names: [id], model_name: model_name))\
            .models.model(model_name).documents.document(record_id)
        end

        def update_records
          request = UpdateRequest.new
          yield request
          magma_client.update(request)
        end
      end

      class FileLinkingWorkflow < Struct.new(:magma_client, :metis_client, :project_name, :bucket_name, :matching_expressions, :attribute_options, keyword_init: true)
        attr_accessor :magma_crud

        def initialize(opts)
          super({attribute_options: {}, matching_expressions: []}.update(opts))
          @magma_crud = MagmaCrudWorkflow.new(magma_client: magma_client)
        end

        def find_matches
          {}.tap do |all_matches|
            metis_client.folders(
                project_name: project_name,
                bucket_name: bucket_name
            ).each do |folder|
              metis_client.list_folder(
                  Etna::Clients::Metis::ListFolderRequest.new(
                      project_name: project_name,
                      bucket_name: bucket_name,
                      folder_path: folder.folder_path,
                  ),
              ).files.each do |file|
                matches = matching_expressions
                          .map { |regex, attribute_name| [regex.match(file.file_path), regex, attribute_name] }
                          .select { |match, regex, attribute_name| !match.nil? }

                if matches.length > 1
                  raise "File path #{file.file_path} matches multiple regex, #{matches.map(&:second)}.  Please modify the matching expressions to disambiguate"
                end

                if matches.length == 1
                  match, _, attribute_name = matches.first
                  match_map = match_map_defaults.dup.update(match.names.zip(match.captures).to_h)
                  key = [match_map, attribute_name]

                  if all_matches.include?(key)
                    raise "Field #{attribute_name} for #{match_map} identified for two files, #{file.file_path} and #{all_matches[key]}.  Please modify the existing matching expressionts to disambiguate"
                  end

                  all_matches[key] = file.file_path
                end
              end
            end
          end
        end

        # Subclasses can override to provide customizable behavior here.
        def match_name_of(model_name)
          model_name
        end

        def match_map_defaults
          { match_name_of("project") => project_name }
        end

        def ensure_containing_record(match_map, model_name)
          raise "Could not find containing model #{model_name} defined." unless (model = models.model(model_name))
          match_name = match_name_of(model_name)

          raise "Matching regex params #{match_map} do not contain #{model_name}'s identifier in key #{match_name}" unless (id = match_map[match_name])
          record = magma_crud.lookup_record(model_name, id)

          if record.nil?
            parent_attribute_name = model.template.attributes.all.select { |a| a.attribute_type == AttributeType::PARENT }.first.attribute_name
            parent_identifier = ensure_containing_record(match_map, model.template.parent)

            attrs = { parent_attribute_name => parent_identifier }
            magma_crud.update_records do |update_request|
              update_request.update_revision(model_name, id, attrs)
            end
          end

          id
        end

        def link_files(model_name)
          find_matches.each do |match_map, attribute_name, file_path|
            id = ensure_containing_record(match_map, model_name)
            magma_crud.update_records do |update_request|
              update_request.update_revision(model_name, id, **revision_for(model_name, id, attribute_name, file_path))
            end
          end
        end

        def revision_for(model_name, id, attribute_name, file_path)
          { attribute_name => { path: "metis://#{project_name}/#{bucket_name}/#{file_path}" } }
        end

        def models
          @models ||= begin
            magma_client.retrieve(RetrievalRequest.new(project_name: self.project_name, model_name: 'all')).models
          end
        end
      end

      # Note!  These workflows are not perfectly atomic, nor perfectly synchronized due to nature of the backend.
      # These primitives are best effort locally synchronized, but cannot defend the backend or simultaneous
      # system updates.
      class ModelSynchronizationWorkflow < Struct.new(:target_client, :source_models, :target_project, keyword_init: true)
        def target_models
          @target_models ||= begin
            target_client.retrieve(RetrievalRequest.new(project_name: self.target_project, model_name: 'all')).models
          end
        end

        def execute_updates(*actions)
          update = UpdateModelRequest.new(project_name: self.target_project)
          actions.each { |a| update.add_action(a) }
          @target_models = nil
          pp update.as_json
          target_client.update_model(update)
        end

        def self.from_api_source(source_project:, source_client:, **kwds)
          self.new(
              source_models: source_client.retrieve(RetrievalRequest.new(project_name: source_project, model_name: 'all')).models,
              **kwds
          )
        end

        # Subclass and override when the source <-> target mapping is not 1:1.
        def target_of_source(model_name)
          model_name
        end

        # Subclass and override when the source <-> target attribute mapping is not 1:1.
        def target_attribute_of_source(model_name, attribute_name)
          attribute_name
        end

        # Potentially cyclical, protected against re-entry via the seen cache.
        # Establishes the link attributes in a given model graph.
        def ensure_model_tree(model_name, seen = Set.new)
          return unless (source_model = source_models.model(model_name))
          return if seen.include?(model_name)
          seen.add(model_name)
          ensure_model(model_name)

          attributes = source_model.template.attributes

          attributes.all.each do |attribute|
            actions = []

            # Don't copy or update parent links.  Once a model has been setup with a parent someway.
            unless attribute.attribute_type == AttributeType::PARENT
              unless attribute.link_model_name.nil?
                ensure_model(attribute.link_model_name)
                next unless (link_model = source_models.model(attribute.link_model_name))
                link_model_attributes = link_model.template.attributes
                reciprocal = link_model_attributes.all.find do |attr|
                  attr.link_model_name == model_name
                end

                actions.push(*add_model_attribute_actions(attribute.link_model_name, reciprocal.attribute_name))
              end

              actions.push(*add_model_attribute_actions(model_name, attribute.attribute_name))
            end
            execute_updates(*actions)

            # Even if it's a parent node, however, we still want to cascade the tree expansion to all links.
            unless attribute.link_model_name.nil?
              ensure_model_tree(attribute.link_model_name, seen)
            end
          end
        end

        def add_model_attribute_actions(model_name, attribute_name)
          return unless (model = source_models.model(model_name))
          return unless (source_attribute = model.template.attributes.attribute(attribute_name))

          target_model_name = target_of_source(model_name)
          target_attribute_name = target_attribute_of_source(model_name, attribute_name)

          target_attributes = target_models.model(target_model_name).template.attributes
          return [] if target_attributes.attribute_keys.include?(attribute_name)

          add_attribute = AddAttributeAction.new(
            model_name: target_model_name,
            attribute_name: target_attribute_name,
          )

          add_attribute.type = source_attribute.attribute_type
          add_attribute.description = source_attribute.desc
          add_attribute.display_name = source_attribute.display_name
          add_attribute.format_hint = source_attribute.format_hint
          add_attribute.hidden = source_attribute.hidden
          add_attribute.link_model_name = source_attribute.link_model_name
          add_attribute.read_only = source_attribute.read_only
          add_attribute.unique = source_attribute.unique
          add_attribute.validation = source_attribute.validation
          add_attribute.restricted = source_attribute.restricted

          [add_attribute]
        end

        # Non cyclical, non re-entrant due to the requirement that parents cannot form a cycle.
        # This method, and it's partner prepare_parent, should never call into any re-entrant or potentially
        # cyclical method, like ensure_model_tree.
        def ensure_model(model_name)
          return unless (source_model = source_models.model(model_name))

          target_model_name = target_of_source(model_name)
          return if target_models.model_keys.include?(target_model_name)

          template = source_model.template

          add_model_action = AddModelAction.new(
              {
                  model_name: target_model_name,
                  identifier: template.identifier,
              }
          )

          parents = template.attributes.all.select { |a| a.attribute_type == AttributeType::PARENT }
          parent_model_name, parent_link_type = prepare_parent(model_name, template, parents)
          unless parent_model_name.nil?
            add_model_action.parent_model_name = parent_model_name
            add_model_action.parent_link_type = parent_link_type
          end

          execute_updates(add_model_action)
        end

        # Non cyclical, non re-entrant due to the requirement that parents cannot form a cycle.
        # This method, and it's partner and ensure_model, should never call into any re-entrant or potentially
        # cyclical method, like ensure_model_tree.
        def prepare_parent(model_name, template, parents)
          return [nil, nil] if parents.empty?
          return [nil, nil] unless (parent_model = source_models.model(template.parent))
          return [nil, nil] unless (child_attribute = parent_model.template.attributes.attribute(model_name))

          ensure_model(template.parent)

          [
              target_of_source(template.parent),
              child_attribute.attribute_type
          ]
        end
      end

      class ShallowCopyModelWorkflow < ModelSynchronizationWorkflow
        def initialize(model_name:, **kwds)
          super(**kwds)
          @model_name = model_name
        end

        # Aside from just creating models for links, do not cascade the expansion.
        def ensure_model_tree(model_name, *args)
          puts "Checking #{model_name} #{@model_name}"
          if model_name == @model_name
            super(model_name, *args)
          end
        end
      end
    end
  end
end