action :create do
  name = new_resource.name
  execute("Create postgres database '#{name}'") do
    user('postgres')
    exists = <<-EOH
      psql --tuples-only -c "select * from pg_database 
where datname='#{name}'" | grep -c #{name}
    EOH
    command("psql -c \"CREATE DATABASE #{name};\"")
    not_if(exists)
  end

  # Lock out the pubic schema
  execute('Lockout postgres schema "public"') do
    user('postgres')
    command("psql -c \"REVOKE ALL ON schema public FROM public;\" -d #{name}")
  end

  # Lock out everyone from the janus DB.
  execute("Lockout public from #{name} DB") do
    user('postgres')
    command("psql -c \"REVOKE ALL ON DATABASE #{name} FROM public;\"")
  end

  # Create the locked down private schema for the janus DB.
  execute("Create \"private\" schema on #{name} DB") do
    user('postgres')
    exists = <<-EOH
      psql -d #{name} -c "SELECT TRUE FROM 
information_schema.schemata WHERE schema_name = 'private';" | grep -c t
    EOH
    command("psql -c \"CREATE SCHEMA private;\" -d #{name}")
    not_if(exists)
  end

  # Allow the 'developer' connection access to the DB.
  execute("Allow 'developer' access to #{name} DB") do
    user('postgres')
    command("psql -c \"GRANT CONNECT ON DATABASE #{name} TO developer;\"")
  end

  # Allow the 'developer' CRUD operations on the DB.
  execute("Allow 'developer' usage on #{name} DB") do
    user('postgres')
    command("psql -d #{name} -c \"GRANT CREATE, USAGE ON SCHEMA private TO developer;\"")
  end
end
