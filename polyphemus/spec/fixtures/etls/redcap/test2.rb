# Define the new classes dynamically
#
# :record, :event, :repeat_instance, :value

define_model("Make").class_eval do
  def identifier(record_name)
    record_name
  end
end

define_model("Model").class_eval do
  def identifier(record_name, event_name)
    event_name
  end
end

define_model("Year").class_eval do
  def identifier(record_name, event_name, (repeat_instrument, repeat_id) )
    [ record_name, event_name, repeat_id ].join(' ')
  end
end
define_model("Feature").class_eval do
  def patch(magma_record_name, record)
    record[:year] = magma_record_name.split('-').values_at(1,2,4).join(' ')
  end
end

def config
  {
    models: {
      make: {
        each: [ "record" ],
        scripts: [
          {
            attributes: {
              date_of_founding: "dof"
            }
          }
        ]
      },
      model: {
        each: [ "record", "event" ],
        scripts: [
          {
            attributes: {
              type: "car_class"
            }
          }
        ]
      },
      year: {
        each: [ "record", "event", "repeat" ],
        scripts: [
          {
            attributes: {
              calendar_year: "year"
            }
          }
        ]
      },
      feature: {
        each: [ "record", "event", "repeat", "value" ],
        scripts: [
          {
            attributes: {
              name: {
                "redcap_field": "feature",
                "value": "select_choice"
              },
              value: "feature"
            }
          }
        ]
      },
    }
  }
end
