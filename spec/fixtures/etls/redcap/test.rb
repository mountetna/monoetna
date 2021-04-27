# Define the new classes dynamically
define_model("ModelOne").class_eval do
  def identifier(record_name, event_name)
    # Hardcode a temp id so that the offset is consistent. Makes
    #   testing less random.
    "::temp-#{record_name}-xyz"
  end
end

define_model("ModelTwo").class_eval do
  def identifier(record_name, event_name)
    record_name
  end
end

define_model("Stats").class_eval do
  def patch(id, record)
    record[:model_two] = id.split("-")[1]
  end
end

def config
  {
    models: {
      model_one: {
        scripts: [
          {
            each: "entity",
            forms: {
              essential_data: {
                birthday: "date_of_birth",
                graduation_date: "commencement",
                name: "name",
              },
            },
          },
        ],
      },
      model_two: {
        scripts: [
          {
            each: "record",
            forms: {
              calendar: {
                label: {
                  redcap_field: "today",
                  value: "label",
                },
                yesterday: {
                  redcap_field: "today",
                  value: "value",
                },
              },
            },
          },
        ],
      },
      stats: {
        scripts: [
          {
            each: "record",
            forms: {
              statistics: {
                height: {
                  redcap_field: "height",
                  value: "value",
                },
                weight: {
                  redcap_field: "weight",
                  value: "value",
                },
              },
            },
          },
        ],
      },
    },
  }
end
