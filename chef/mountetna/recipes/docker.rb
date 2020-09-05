# docker_service 'default' do
#   # package_options %q|--nobest| # Issue with runc in containerd right now for RHEL8
# end

docker = node['docker']
swarm = docker['swarm']

addr = swarm['leader_addr']
token = swarm['token']
registry = docker['registry']

url = registry['url']
username = registry['username']
password = registry['password']
email = registry['email']

docker_installation_package 'default' do
  version '19.03.12'
  action :create
end

docker_service_manager 'default' do
  action :start
  storage_driver 'overlay2'
  data_root docker['data_root']
  host docker['host']
end

if node['is_kitchen']
  # HACK!  We can't fully run docker in docker without a lot of other work,
  # So we actually inject the host docker's socket into the container, but
  # it seems that installing docker is 'too helpful' and clears any existing
  # docker.socket it finds from previous installations.  This hack, used by kitchen
  # links in the injected docker.sock AFTER installation
  execute 'is_kitchen hack to link docker.sock' do
    user("root")
    command("ln -s /usr/var/run/docker.sock /run/docker.sock")
    not_if { ::File.exists?('/run/docker.sock')}
  end
end

if [username, password, email].all? { |s| s && !s.empty? }
  docker_registry url do
    username username
    password password
    email email
  end
else
  Chef::Log.warn("docker username, password, and email were not set, cannot complete registry login.")
end

if [token, addr].all? { |s| s && !s.empty? }
  execute "Join swarm at #{addr}" do
    user("root")
    command("docker swarm join --token #{token} #{addr}")
  end
else
  Chef::Log.warn("leader_addr and token required to complete swarm join, this node will not be available.")
end

mountetna_systemd_wrapped_container 'watchtower' do
  def watchtower_notification_options
    slack = node['slack_notifications']
    watchtower = slack['watchtower']
    hook_url = watchtower['hook_url']
    channel = watchtower['channel']

    return [] unless hook_url && !hook_url.empty?

    result = []
    hostname=`hostname`.strip
    result << "-e 'WATCHTOWER_NOTIFICATIONS=slack'"
    result << "-e 'WATCHTOWER_NOTIFICATION_SLACK_HOOK_URL=#{hook_url}'"
    result << "-e 'WATCHTOWER_NOTIFICATION_SLACK_IDENTIFIER=Watchtower running on #{hostname}'"
    result << "-e 'WATCHTOWER_NOTIFICATION_SLACK_CHANNEL=#{channel}'"
    result << "-e 'WATCHTOWER_NOTIFICATION_SLACK_ICON_EMOJI=:whale:'"
  end

  image 'containrrr/watchtower'
  options [
      "-v /var/run/docker.sock:/var/run/docker.sock",
  ] + watchtower_notification_options
end

