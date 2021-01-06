# Define the new classes dynamically
define_model("ModelOne")
define_model("ModelTwo")

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
