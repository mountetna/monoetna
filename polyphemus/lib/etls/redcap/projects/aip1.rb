# Define the new classes dynamically

# In order to get consistent date-shifting, we'll include
#   this module.
module SubjectTable
  def self.included(base)
    base.class_eval do
      # We override these methods via the included hook,
      #   otherwise just including the SubjectTable module
      #   means they will be overwritten by the default
      #   method definitions.
      def identifier(record_name, event_name = nil, identifier_fields: nil)
        raise ArgumentError, "Missing autoipi_id field in model\'s identifier_fields config" unless !identifier_fields.nil? && !identifier_fields[:autoipi_id].nil?

        [
          "::temp", identifier_fields[:autoipi_id], rand(36 ** 8).to_s(36),
        ].compact.join("-")
      end

      def subject_id(magma_temp_id)
        # because the autoipi_ids have `-` in them, we'll slice out
        #   the relevant parts and re-join...
        magma_temp_id.split(/-/)[1..-2].join("-")
      end

      def offset_id(magma_record_name)
        subject_id(magma_record_name)
      end

      def patch(magma_record_name, record)
        set_record_subject(magma_record_name, record)
      end

      def set_record_subject(magma_record_name, record)
        record[:subject] = subject_id(magma_record_name)
      end
    end
  end
end

define_model("Subject").class_eval do
  def identifier(record_name, event_name = nil, identifier_fields: nil)
    raise ArgumentError, "Missing autoipi_id field from form" unless !identifier_fields.nil? && !identifier_fields[:autoipi_id].nil?

    identifier_fields[:autoipi_id]
  end

  def patch(magma_record_name, record)
    record[:current_age] = [89, record[:current_age]].min if record[:current_age]
  end

  def offset_id(magma_record_name)
    magma_record_name
  end
end

define_model("AutoimmuneHistory").class_eval do
  include SubjectTable

  def patch(id, record)
    super

    record[:age] = [89, record[:age]].min if record[:age]
  end
end

define_model("Medication").class_eval do
  include SubjectTable
end

define_model("Comorbidity").class_eval do
  include SubjectTable
end

define_model("Demographic").class_eval do
  include SubjectTable
end
