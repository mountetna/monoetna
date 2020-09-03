# docker_service 'default' do
#   # package_options %q|--nobest| # Issue with runc in containerd right now for RHEL8
# end

docker_installation_package 'default' do
  version '19.03.12'
  action :create
end

docker_service_manager 'default' do
  action :start
  storage_driver 'overlay2'
  data_root node['docker']['data_root']
end
