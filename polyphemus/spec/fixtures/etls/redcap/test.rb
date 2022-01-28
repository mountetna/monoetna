# Define the new classes dynamically
define_model("ModelOne").class_eval do
  def identifier(record_name, identifier_fields: nil)
    # Hardcode a temp id so that the offset is consistent. Makes
    #   testing less random.
    "::temp-#{record_name}-xyz"
  end

  def offset_id(record_name)
    record_name
  end
end

define_model("ModelTwo").class_eval do
  def identifier(record_name, identifier_fields: nil)
    record_name
  end

  def offset_id(record_name)
    record_name
  end

  def patch(id, record)
    record[:model_one] = id + "-one"
  end

  def ensure_containing_records?
    true
  end

  def ensure_records(current_records)
    model_one_names = current_records[:model_two].values.map { |r| r[:model_one] }.uniq.compact

    result = current_records[:model_one].dup || {}
    
    model_one_names.map do |model_one_name|
      result[model_one_name] ||= {}
      result[model_one_name][:parent_model] = "#{model_one_name}-parent"
    end

    {
      model_one: result
    }
  end
end

define_model("Stats").class_eval do
  def patch(id, record)
    record[:model_two] = id.split("-")[1]
  end

  def offset_id(record_name)
    record_name
  end
end

define_model("Citation").class_eval do
  def redcap_id(record_name, record)
    record_name.split("-")[1..2]
  end

  def offset_id(record_name)
    record_name
  end
end

define_model("BadModel").class_eval do
  def identifier(record_name, identifier_fields: nil)
    # Hardcode a temp id so that the offset is consistent. Makes
    #   testing less random.
    "::temp-#{record_name}-abc"
  end
end

define_model("ModelWithAlternateId").class_eval do
  def identifier(record_name, identifier_fields: nil)
    # Hardcode a temp id so that the offset is consistent. Makes
    #   testing less random.
    raise ArgumentError, "Missing :date_of_birth in form" if identifier_fields.nil? || identifier_fields[:date_of_birth].nil?
  
    "::temp-#{identifier_fields[:date_of_birth]}-abc"
  end
end
