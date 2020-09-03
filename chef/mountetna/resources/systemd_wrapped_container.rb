property :image, String, required: true
property :tag, String, default: 'latest'
property :cmd, String, default: ''
property :options, Array, default: []
property :Unit, Hash, default: {}
property :Service, Hash, default: {}
property :start, [TrueClass, FalseClass], default: true
property :restart, default: 'always'

action :create do
  docker_image new_resource.image do
    action :pull
    tag new_resource.tag
  end

  systemd_unit "#{new_resource.name}.service" do
    content({
        Unit: {
            Description: "Runs #{new_resource.name} inside docker image #{new_resource.image}",
            After: "docker.service",
            Requires: "docker.service",
        }.update(new_resource.Unit),
        Service: {
            Type: 'simple',
            User: 'root',
            WorkingDirectory: '/root',
            Environment: "",
            Restart: new_resource.restart,
            ExecStart: "/usr/bin/docker run --name #{new_resource.name} --rm #{new_resource.options.join(' ')} #{new_resource.image}:#{new_resource.tag} #{new_resource.cmd}",
        }.update(new_resource.Service),
    })

    action([:create, :enable] + (new_resource.start ? [:start] : []))
  end
end
