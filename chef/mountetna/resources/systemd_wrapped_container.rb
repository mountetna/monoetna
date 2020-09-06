property :image, String, required: true
property :tag, String, default: 'latest'
property :cmd, String, default: ''
property :options, Array, default: []
property :Unit, Hash, default: {}
property :Service, Hash, default: {}
property :start, [TrueClass, FalseClass], default: true
property :restart, default: 'always'

action :create do
  container_resource = new_resource

  docker_image container_resource.image do
    action :pull
    tag container_resource.tag
  end

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
            Restart: container_resource.restart,
            ExecStart: "/usr/bin/docker run --name #{container_resource.name} --rm #{container_resource.options.join(' ')} #{container_resource.image}:#{container_resource.tag} #{container_resource.cmd}",
        }.update(container_resource.Service),
    })

    action([:create, :enable] + (container_resource.start ? [:start] : []))
  end
end
