# Define the new classes dynamically
define_model("ModelOne")
define_model("ModelTwo").class_eval do
  def identifier(record_name, event_name)
    record_name
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
                name: "name"
              }
            }
          }
        ]
      },
      model_two: {
        scripts: [
          {
            each: "record",
            forms: {
              calendar: {
                label: {
                  redcap_field: "today",
                  value: "label"
                },
                yesterday: {
                  redcap_field: "today",
                  value: "value",
                  exists: true
                }
              }
            }
          }
        ]
      }
    }
  }
end
