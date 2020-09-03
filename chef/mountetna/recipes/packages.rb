# Installs global system packages.
# NOTE: Prefer keeping application related packages inside the related application containers.
# System level packages should be small and operationally useful.
# These packages will be included in all instances that run mountetna, so they should be valuable and necessary.
[
    'git',
    'tmux',
    'openssl',
    'curl',
    'postgresql-devel',
    'postgresql-contrib',
    'postgresql-server',
    'vim',
].each do |p|
  package(p) do
    action(:install)
  end
end
