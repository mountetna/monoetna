property :enabled, [ TrueClass, FalseClass ], default: false
property :has_db, [ TrueClass, FalseClass ], default: true
property :docker_tag, String, required: true

action :create do
  base_v2_app_db name if has_db
end
