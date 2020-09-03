def postgres_initialized?
  base = '/var/lib/pgsql/data/base'
  global = '/var/lib/pgsql/data/global'
  File.directory?(base) && File.directory?(global)
end

unless postgres_initialized?
  execute('Initialize a postgres cluster') do
    user('root')
    command('postgresql-setup initdb')
  end
end

template('/var/lib/pgsql/data/pg_hba.conf') do
  source('pg_hba.conf.erb')
  owner('postgres')
  group('postgres')
  mode('600')
end

execute('Start postgres') do
  user('root')
  command('systemctl start postgresql')
end

execute 'Enable postgres' do
  user('root')
  command('systemctl enable postgresql')
end

execute('Create postgres role "developer"') do
  user('postgres')
  exists = <<-EOH
      psql -c "SELECT 1 FROM pg_roles WHERE 
rolname=\'developer\'" | grep -q 1
  EOH
  command('psql -c "CREATE ROLE developer WITH LOGIN;"')
  not_if(exists)
end

if node['psql_developer_password'].nil?
  fail "node['psql_developer_password'] is nil, but requires a value."
end

execute('Lock the "developer" postgres user') do
  user('postgres')
  psql_developer_password = node['psql_developer_password']
  command(%Q(psql -c "ALTER USER developer WITH ENCRYPTED PASSWORD \'#{psql_developer_password}\'";))
end
