property :has_db, [TrueClass, FalseClass], default: true
property :tag, String, default: node['docker']['default_tag']
property :extra_options, Array, default: []

action :create do
  name = new_resource.name

  mountetna_systemd_wrapped_container "#{name}_app" do
    image "mountetna/#{name}"
    tag new_resource.tag
    cmd "./bin/puma.sh"

    options [
        "-e 'APP_NAME=#{name}'",
        "-e '#{name.upcase}_ENV=production",
    ] + new_resource.extra_options
  end

  mountetna_systemd_wrapped_container "#{name}_app_fe" do
    image "mountetna/etna-apache"
    tag new_resource.tag
    cmd "./bin/puma.sh"

    options [
        "-e 'APP_NAME=#{name}'",
        "-e '#{name.upcase}_ENV=production",
        "--link #{name}_app:#{name}_app", # The link tells watchtower to bring this service down and restart it after pulling updates to the original
        "--volumes-from #{name}_app"
    ]
  end

  mountetna_app_db name if new_resource.has_db
end