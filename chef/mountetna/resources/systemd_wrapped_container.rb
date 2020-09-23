property :image, String, required: true
property :tag, String, default: 'latest'
property :cmd, String, default: ''
property :options, Array, default: []
property :Unit, Hash, default: {}
property :Service, Hash, default: {}
property :start, [TrueClass, FalseClass], default: true
property :extra_setup, Array

action :create do
  container_resource = new_resource

  docker_image container_resource.image do
    action :pull
    tag container_resource.tag
  end

  container_options = "--name #{container_resource.name} #{container_resource.options.join(' ')} #{container_resource.image}:#{container_resource.tag} #{container_resource.cmd}"

  # execute "stopping existing container #{container_resource.name}" do
  #   command "docker stop #{container_resource.name}"
  #   only_if { `docker ps`.include? "#{container_resource.name}_app" }
  # end
  #
  # execute "removing existing container #{container_resource.name}" do
  #   command "docker rm #{container_resource.name}"
  #   only_if { `docker ps -all`.include? "#{container_resource.name}_app" }
  # end
  #
  # execute "create docker container #{container_resource.name}" do
  #   command "docker container create #{container_options}"
  #   not_if { `docker ps --all`.include? "#{container_resource.name}_app" }
  # end

  systemd_unit "#{container_resource.name}.service" do
    content({
        Unit: {
            Description: "Runs #{container_resource.name} inside docker image #{container_resource.image}",
            After: "docker.service",
            Requires: "docker.service",
        }.update(container_resource.Unit),
        Service: {
            Type: 'simple',
            User: 'root',
            WorkingDirectory: '/root',
            Environment: "",
            Restart: 'always',
            RestartSec: '3',
            ExecStart: "/usr/bin/docker run --rm #{container_options}",
        }.update(container_resource.Service),
    })

    action([:create, :enable] + (container_resource.start ? [:restart] : []))
  end

  unless container_resource.extra_setup.nil?
    container_resource.extra_setup.each do |c|
      execute "running #{c} for #{container_resource.name}" do
        command c
      end
    end
  end
end
