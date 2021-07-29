# Define the new classes dynamically
define_model("ModelOne").class_eval do
  def identifier(record_name)
    # Hardcode a temp id so that the offset is consistent. Makes
    #   testing less random.
    "::temp-#{record_name}-xyz"
  end
end

define_model("ModelTwo").class_eval do
  def identifier(record_name)
    record_name
  end
end

define_model("Stats").class_eval do
  def patch(id, record)
    record[:model_two] = id.split("-")[1]
  end
end

define_model("Citation").class_eval do
  def redcap_id(record_name, record)
    record_name.split('-')[1..2]
  end
end

def config
  {
    models: {
      model_one: {
        each: [ "record" ],
        scripts: [
          {
            attributes: {
              birthday: "date_of_birth",
              graduation_date: "commencement",
              name: "name"
            }
          }
        ]
      },
      model_two: {
        each: [ "record" ],
        scripts: [
          {
            attributes: {
              label: {
                redcap_field: "today",
                value: "label"
              },
              yesterday: {
                redcap_field: "today",
                value: "value"
              }
            }
          }
        ]
      },
      stats: {
        each: [ "record" ],
        scripts: [
          {
            attributes: {
              height: {
                redcap_field: "height",
                value: "value"
              },
              weight: {
                redcap_field: "weight",
                value: "value"
              }
            }
          }
        ]
      },
      citation: {
        each: [ "record", "event" ],
        invert: true,
        scripts: [
          {
            attributes: {
              name: {
                match: /citation-(111|abc)/,
                value: "none"
              },
              date: {
                redcap_field: "citation_date",
                value: "value"
              }
            }
          }
        ]
      }
    }
  }
end
