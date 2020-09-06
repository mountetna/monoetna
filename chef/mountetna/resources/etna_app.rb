require 'yaml'

property :has_db, [TrueClass, FalseClass], default: true
property :tag, String, default: node['docker']['default_tag']
property :extra_docker_options, Array, default: []
property :extra_yml_config, Hash, default: {}

action :create do
  name = new_resource.name

  directory "/var/mountetna/#{name}" do
    recursive true
    owner "root"
    group 'docker'
  end

  file "/var/mountetna/#{name}/config.yml" do
    owner 'root'
    group 'docker'
    content(EtnaConfigBuilder.new(name, node).merge(new_resource.extra_yml_config).build.to_yaml)
    # Can contain secrets; ensure only true administrative users can read this file.
    mode '640'
  end

  template "/var/mountetna/httpd.conf.d/#{name}.app.include" do
    owner 'root'
    group 'docker'
    source 'etna_app.include.erb'

    variables(
        app_name: name,
        host_prefix: URI.parse(node['hosts'][name]).host.gsub('.', '_')
    )
    mode '644'
  end

  mountetna_systemd_wrapped_container "#{name}_app" do
    image "etnaagent/#{name}"
    tag new_resource.tag
    cmd "./bin/puma.sh"

    options [
        "-e 'APP_NAME=#{name}'",
        "-e '#{name.upcase}_ENV=production'",
        "--mount 'source=/var/mountetna/#{name}/config.yml,target=/app/config.yml,readonly'"
    ] + new_resource.extra_docker_options
  end

  mountetna_systemd_wrapped_container "#{name}_app_fe" do
    image "etnaagent/etna-apache"
    tag new_resource.tag
    cmd "./bin/puma.sh"

    options [
        "-e 'APP_NAME=#{name}'",
        "-e '#{name.upcase}_ENV=production'",
        "--link #{name}_app:#{name}_app", # The link tells watchtower to bring this service down and restart it after pulling updates to the original
        "--volumes-from #{name}_app"
    ]
  end

  mountetna_app_db(name) if new_resource.has_db
end

class EtnaConfigBuilder
  def initialize(name, node)
    config = @config = {}
    production_config = config[:production] = {}
    production_config[:auth_redirect] = node['hosts']['janus']
    production_config[:token_algo] = node['janus_token']['algo']
    production_config[:token_name] = node['janus_token']['name']
    # With mountetna, logs will redirect to the journalctl process running them.
    production_config[:log_file] = '/dev/stdout'
    production_config[:log_level] = 'info'
    production_config[:rsa_public] = node['janus_token']['rsa_public']

    # Add in all the hosts.
    node['hosts'].each do |k, v|
      (production_config[k.to_sym] ||= {})[:host] = v
    end

    db_config = production_config[:db] = {}
    db_config[:host] = '192.168.0.1'

    db_config[:adapter] = 'postgesql'
    db_config[:encoding] = 'unicode'
    db_config[:user] = 'developer'
    db_config[:password] = node['psql_developer_password']
    db_config[:pool] = 5
    db_config[:timeout] = 5000
    db_config[:database] = name

    rollbar_config = (production_config[:rollbar] ||= {})
    rollbar_config[:access_token] = node['rollbar']['access_token']
  end

  def merge(extra)
    extra.each do |k, v|
      if v.is_a? Hash
        if @config.include?(k) && !@config[k].is_a?(Hash)
          raise "Key #{k} in extra_yml_config conflicts with etna_app's base config."
        end
        (@config[k] ||= {}).update(v)
      else
        @config[k] = v
      end
    end

    self
  end

  def ensure_no_nils(yml, path = [])
    yml.each do |k, v|
      if v.is_a? Hash
        ensure_no_nils(v, path + [k])
      end

      raise "Unexpected nil in etna_app config found: #{path.join('.')}.#{k}" if v.nil?
    end

    yml
  end

  def build
    ensure_no_nils(@config)
  end
end