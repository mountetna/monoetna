require_relative './loader/temp_id'
require_relative './loader/multi_update'
require_relative './loader/record_entry'
require_relative './loader/record_hierarchy_cache'

class Magma
  class LoadFailed < Exception
    attr_reader :complaints

    def initialize(complaints)
      @complaints = complaints
    end
  end

  class BaseAttributeEntry
    def initialize(model, attribute, loader)
      @model = model
      @attribute = attribute
      @loader = loader
    end

    def entry value
      nil
    end
  end

  # A generic loader class.
  class Loader
    class << self
      def description desc=nil
        @description ||= desc
      end

      def loader_name
        name.snake_case.sub(/_loader$/,'')
      end
    end

    attr_reader :validator, :user, :dry_run

    def initialize(user, project_name, dry_run: false)
      @user = user
      @project_name = project_name
      @validator = Magma::Validation.new
      @censor = Magma::Censor.new(@user,@project_name)
      @records = {}
      @temp_id_counter = 0
      @attribute_entries = {}
      @identifiers = {}
      @now = Time.now.iso8601
      @dry_run = dry_run
      @flags = Magma::Project.flags(project_name)
    end

    def push_record(model, record_name, revision)
      records(model)[record_name] ||= RecordEntry.new(model, record_name, self)

      records(model)[record_name] << revision
    end

    def to_payload
      payload = Magma::Payload.new

      @records.each do |model, record_set|
        next if record_set.empty?

        payload.add_model(model)

        payload.add_records(model, record_set.values.map(&:payload_entry))
      end

      return payload.to_hash
    end

    def is_collection_attribute?(attribute)
      attribute.is_a?(Magma::CollectionAttribute)
    end

    def is_child_attribute?(attribute)
      attribute.is_a?(Magma::ChildAttribute)
    end

    def is_table_attribute?(attribute)
      attribute.is_a?(Magma::TableAttribute)
    end

    def is_incoming_link?(attribute)
      # Based on the attribute class, determines if this is the incoming attribute,
      #   for a foreign-key relationship. If so, return true.
      is_collection_attribute?(attribute) ||
      is_child_attribute?(attribute) ||
      is_table_attribute?(attribute)
    end

    def explicit_child_revision_exists?(revisions, model_name, record_name, attribute_name)
      # Check if the child's parent / link is being explicitly set by the user
      !!revisions.dig(model_name.to_sym, record_name.to_sym, attribute_name.to_sym)
    end

    def explicit_link_revision_exists?(revisions, model_name, record_name, attribute, child_record_name)
      # We need to check if the child record has
      #   been set for ANY record in the same parent_model.

      # NOTE: different check behavior for CollectionAttribute vs. ChildAttribute!
      explicit_revisions_found = 0

      revisions[model_name.to_sym].each do |record_name, record_revisions|
        next unless record_revisions.key?(attribute.attribute_name.to_sym)

        attribute_revision = record_revisions[attribute.attribute_name.to_sym]

        explicit_revisions_found += 1 if is_collection_attribute?(attribute) && attribute_revision.include?(child_record_name)
        explicit_revisions_found += 1 if is_child_attribute?(attribute) && child_record_name == attribute_revision
      end

      explicit_revisions_found > 0
    end

    def explicit_revision_exists?(revisions:, parent_model:, parent_record_name:, child_record_name:, parent_attribute:)
      # Cannot find explicit revisions from @records, because
      #   some of those are calculated! So we have to look for
      #   explicit revisions from the user-supplied revisions hash.
      # Note that because explicit revisions can be either parent / link -> child
      #   or child -> parent / link, we need to check both cases.
      explicit_child_revision_exists?(
        revisions,
        parent_attribute.link_model.model_name,
        child_record_name,
        parent_attribute.link_attribute_name) ||
      explicit_link_revision_exists?(
        revisions,
        parent_model.model_name,
        parent_record_name,
        parent_attribute,
        child_record_name)
    end

    def find_parent_models(identifier)
      # This method attempts to find all parent models and (their identifiers) for a given identifier.
      #
      # It returns an ORDERED list of hashes, where each hash contains the model, the identifier, and the
      # parent_model. The order corresponds to the model hierarchy - and the project model is the first item in the list.
      # If a identifier does not exist in the hierarchy chain, an empty list is returned.
      #
      # This method succeeds at finding parent identifiers if three conditions are met:
      # 1. A proper grammar is created with gnomon and the grammar rules correspond to a model hierarchy.
      # 2. All identifiers in the hierarchy have been created in the identifier table.
      # 3. All models have been created in the db.
      #
      grammar = Magma::Gnomon::Grammar.for_project(@project_name)
      grammar_decomposed = grammar.decompose(identifier)
      return [] if grammar_decomposed.nil?

      # First we gather all available parent models and identifiers associated with the identifier.
      # If someone has created a well formatted grammar, with proper identifiers,
      # there MAY be a model hierarchy chain we can infer
      available_models = []
      grammar_decomposed[:rules].each do |token, token_hash|
        begin
          model = Magma.instance.get_model(@project_name, token)
          available_models << {"model": token.to_sym, "identifier": token_hash[:name], "parent_model": model.parent_model_name}
        rescue NameError => e
          next if e
        end
      end

      # Next we determine if the hierarchy chain is valid for the models that we've found.
      partitioned = available_models.partition { |h| h[:parent_model].nil? }
      root_model, other_models = partitioned[0][0], partitioned[1].flatten
      hierarchy_count = other_models.count

      # We do this by recursively searching for children top-down.
      # We start with the projects model, and then search for any models that have the key
      # :parent_model = :projects. Lets say we find some model X, then the following iteration,
      # this repeats, and we search for any models where :parent_model = :X. Each model that is found
      # is appended to the sorted_array.
      def sort(arr, sorted_array, model_name)
        # Look for the child model
        child_model = arr.find { |hash| hash[:parent_model] == model_name }
        if child_model and !arr.empty?
          sorted_array << child_model
          arr.delete(child_model)
          sort(arr, sorted_array, child_model[:model])
        else
          return
        end
      end

      sorted = []
      sort(other_models, sorted, root_model[:model])

      # If we successfully infer a model hierarchy, the counts should be the same
      hierarchy_count == sorted.count ? sorted.prepend(root_model) : []
    end

    def push_parent_identifiers
      # If a project is configured with the 'identifier' flag, anytime a record with an identifier is created or updated,
      # we can implicitly infer what their parents must be via find_parent_models(). We then update the record and its revision
      # to include its foreign key, and do this up the model hierarchy.
      #
      # This method handles a very specific update:
      # - The incoming update does NOT explicitly reference a parent/child's foreign key.
      # - It doesn't matter what type of attributes we are updating
      @records.each do |model, record_set|
        record_set.each do |record_name, record|
          gnomon_mode = Magma::Flags::GNOMON_MODE
          flag_value = @flags[gnomon_mode[:name]]

          next unless flag_value == gnomon_mode[:identifier] and not record.includes_parent_record?

          # Attempt to find parent models
          parent_models = find_parent_models(record_name)

          next if parent_models.empty?

          # Update foreign keys
          parent_models.each do |parent|
            # TODO: push_records() starting from parents
            # TODO: update foreign key of original record
          end

        end
      end
    end

    def push_implicit_link_revisions(revisions)
      # When updating link or parent attributes from the top-down,
      #   we may not know what the previous relationships were.
      # For example, changing a link child from record A to B,
      #   the explicit revision is LinkModel -> B.
      # But there is also an implicit revision, of updating
      #   record A to have a `nil` parent.
      # So we also have to push records for all the
      #   implicit link revisions, when the attribute
      #   is a Child or Collection type (i.e. the revision
      #   comes from the parent / link).
      # But, in a multi-revision scenario, we should
      #   only push the implicit revisions when those records + attributes
      #   themselves aren't being revised, otherwise we risk
      #   overwriting an explicit revision.

      # We iterate over revisions to see what all has been updated.
      # If there are any ChildAttribute or CollectionAttribute values
      #   that changed, we'll need to investigate further if any implicit
      #   revisions exist.
      # Do not do this over @records, because some of those revisions
      #   are calculated and could lead to incorrectly disconnecting
      #   currently-attached records.

      revisions.each do |model_name, model_revisions|
        model = Magma.instance.get_model(@project_name, model_name)

        # For each model, collect all the record_names being revised
        #   for each child / collection attribute.
        revised_model_records = {}

        model_revisions.each do |record_name, revision|
          revision.each do |attribute_name, value|
            attribute = model.attributes[attribute_name]

            if is_collection_attribute?(attribute) || is_child_attribute?(attribute)
              revised_model_records[attribute] ||= []
              revised_model_records[attribute] << record_name.to_s
            end
          end
        end
        implicit_revisions(revisions, model, revised_model_records) do |child_model, record_name, attribute_name|
          push_record(
            child_model, record_name.to_s,
            attribute_name.to_sym => nil,
            updated_at: @now)
        end
      end
    end

    def implicit_revisions(revisions, parent_model, revised_model_records)
      # Here we fetch all current records that have one of the parent model::record_names as the
      #   parent, and exclude any the user is explicitly setting in `revisions`.
      revised_model_records.keys.each do |attribute|
        child_model = attribute.link_model
        parent_record_names = revised_model_records[attribute]

        question = Magma::Question.new(@project_name, [
          child_model.model_name,
          [attribute.link_attribute_name, '::identifier', '::in', parent_record_names],
            '::all', attribute.link_attribute_name, '::identifier'
        ], user: @user)
        current_record_names = question.answer.to_a

        current_record_names.reject { |child_record_name, parent_record_name|
          explicit_revision_exists?(
            revisions: revisions,
            parent_model: parent_model,
            parent_record_name: parent_record_name,
            child_record_name: child_record_name,
            parent_attribute: attribute)
         }.each do |child_record_name, parent_record_name|
          yield [child_model, child_record_name, attribute.link_attribute_name]
        end
      end
    end

    def push_links(model, record_name, revision)
      revision.each do |attribute_name, value|
        next unless model.has_attribute?(attribute_name)

        attribute = model.attributes[attribute_name]

        attribute.revision_to_links(record_name, value) do |link_model, link_identifiers|
          # When the explicit revision is from the parent / link -> children, the
          #   new link records need to be single-entry records,
          #   because they are linking from child up to a single
          #   parent.
          # When the explicit revision is from the child -> parent / link,
          #   the new link records need to be Arrays,
          #   because they are linking from parent to a collection
          #   of records.
          link_record_name = is_incoming_link?(attribute) ?
            record_name.to_s :
            [ record_name.to_s ]

          link_identifiers&.each do |link_identifier|
            next if link_identifier.nil?
            push_record(
              link_model, link_identifier,
              attribute.link_attribute_name.to_sym => link_record_name,
              created_at: @now,
              updated_at: @now
            )
          end
        end
      end
    end

    # Once we have loaded up all the records we wish to insert/update (upsert)
    # we run this function to kick off the DB insert and update queries.
    def dispatch_record_set
      validate!

      censor_revisions!

      run_attribute_hooks!

      begin
        upsert

        update_temp_ids
      end unless @dry_run

      payload = to_payload

      reset

      return payload
    end

    def reset
      @records = {}
      @validator = Magma::Validation.new
      @censor = Magma::Censor.new(@user,@project_name)
      @attribute_entries = {}
      @identifiers = {}
      GC.start
    end

    # This lets you give an arbitrary object (e.g. a model used in the loader) a
    # temporary id so you can make database associations.
    def temp_id(obj)
      return nil if obj.nil?
      temp_ids[obj] ||= TempId.new(new_temp_id, obj)
    end

    TEMP_ID_MATCH=/^::/

    def identifier_id(model, identifier)
      return nil unless identifier

      @identifiers[model] ||= model.select_map(
        [model.identity.column_name.to_sym, :id]
      ).map do |identifier, id|
        [ identifier.to_s, id ]
      end.to_h

      @identifiers[model][identifier] ||= temp_id([ model, identifier ])

      @identifiers[model][identifier]
    end

    def real_id(model, identifier)
      return identifier if @dry_run

      id = identifier_id(model, identifier)
      id.is_a?(TempId) ? id.real_id : id
    end

    def identifier_exists?(model, identifier)
      !identifier_id(model, identifier).is_a?(TempId)
    end

    def records(model)
      return @records[model] if @records[model]

      @records[model] = {}

      ensure_link_models(model)

      @records[model]
    end

    # This is implemented in the loader because it
    #   requires access to both the set of @records as well as the database.
    def path_to_date_shift_root(model, record_name)
      @record_hierarchy_cache ||= Magma::RecordHierarchyCache.new(@records)

      @record_hierarchy_cache.path_to_date_shift_root(model, record_name)
    end

    def is_connected_to_date_shift_root?(model, record_name)
      !path_to_date_shift_root(model, record_name).empty?
    end

    private

    def validate!
      complaints = []

      @records.each do |model, record_set|
        next if record_set.empty?
        complaints.concat(record_set.values.map(&:complaints))
      end

      complaints.flatten!

      raise Magma::LoadFailed.new(complaints) unless complaints.empty?
    end

    def censor_revisions!
      complaints = []
      @records.each do |model, record_set|
        reasons = @censor.censored_reasons(model, record_set)

        next if reasons.empty?

        complaints += reasons
      end

      raise Magma::LoadFailed.new(complaints) unless complaints.empty?
    end

    def run_attribute_hooks!
      bulk_load_type = {}
      @records.each do |model, record_set|
        model.attributes.each do |att_name, attribute|
          bulk_load_attribute = {}

          record_set.each do |record_name, record|
            next unless record.has_key?(att_name)
            error = attribute.load_hook(self, record_name, record[att_name], bulk_load_attribute)

            raise Magma::LoadFailed.new([error]) if error

            error = attribute.patch_record_load_hook(self, record_name, record)

            raise Magma::LoadFailed.new([error]) if error
          end

          error = attribute.bulk_load_hook(self, bulk_load_attribute)
          raise Magma::LoadFailed.new([error]) if error

          unless bulk_load_attribute.empty?
            bulk_load_type[ attribute.class ] ||= {}
            bulk_load_type[ attribute.class ][ attribute ] = bulk_load_attribute
          end
        end
      end

      bulk_load_type.each do |attribute_class, bulk_type_attributes|
        error = attribute_class.type_bulk_load_hook(self, @project_name, bulk_type_attributes)
        raise Magma::LoadFailed.new([error]) if error
      end
    end

    # This 'upsert' function will look at the records and either insert or
    # update them as necessary.
    def upsert
      # Loop the records separate them into an insert group and an update group.
      # @records is separated out by model.
      @records.each do |model, record_set|
        # Skip if the record_set for this model is empty.
        next if record_set.empty?

        # Our insert and update record groupings.
        insert_records = record_set.values.select(&:valid_new_entry?)
        update_records = record_set.values.select(&:valid_update_entry?)

        # Run the record insertion.
        multi_insert(model, insert_records)

        # Run the record updates.
        multi_update(model, update_records)
      end
    rescue Exception => e
      raise Magma::LoadFailed.new([e.message])
    end

    def multi_insert(model, insert_records)
      by_attribute_key(insert_records) do |records|
        insert_ids = model.multi_insert(
          records.map(&:insert_entry),
          return: :primary_key
        )

        if insert_ids
          records.zip(insert_ids).each do |record, real_id|
            record.real_id = real_id
          end
        end
      end
    end

    def multi_update(model, update_records)
      by_attribute_key(update_records) do |records|
        MultiUpdate.new(model, records.map(&:update_entry), :id, :id).update
      end
    end

    def by_attribute_key(all_records)
      all_records.group_by(&:attribute_key).each do |_, records|
        yield records
      end
    end

    def update_temp_ids
      @records.each do |model, record_set|
        next if record_set.empty?
        temp_records = record_set.values.select(&:valid_temp_update?)

        MultiUpdate.new(model, temp_records.map(&:temp_entry), :real_id, :id).update
      end
    end

    def temp_ids
      @temp_ids ||= {}
    end

    def new_temp_id
      @temp_id_counter += 1
    end

    def ensure_link_models(model)
      model.attributes.each do |att_name, att|
        records(att.link_model) if att.is_a?(Magma::Link)
      end
    end
  end
end
