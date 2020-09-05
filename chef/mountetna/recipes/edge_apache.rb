include_recipe 'mountetna::docker'

directory '/var/mountetna/httpd.conf.d' do
  recursive true
  owner "root"
  group 'docker'
end

cookbook_file '/var/mountetna/httpd.conf.d/main.conf' do
  owner 'root'
  group 'docker'
  source 'httpd.conf.d/main.conf'
  mode '640'
end


# TODO: Need to mount /app/data directory to serve metis assets from.  When we get there.
mountetna_systemd_wrapped_container 'edge-apache' do
  image 'etnaagent/edge-apache'
  tag node['docker']['default_tag']
  options [
      "-v /var/mountetna/edge-apache/httpd.conf.d/:/usr/opt/httpd.conf.d/:ro"
  ]
end

