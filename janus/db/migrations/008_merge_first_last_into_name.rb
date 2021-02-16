Sequel.migration do 
  up do
    alter_table(:users) do
      add_column :name, String
    end

    Janus.instance.db.execute("UPDATE users SET name=CONCAT(first_name,' ', last_name);")

    alter_table(:users) do
      set_column_not_null :name
      drop_column :first_name
      drop_column :last_name
    end
  end

  down do
    alter_table(:users) do
      add_column :first_name, String
      add_column :last_name, String
    end

    Janus.instance.db.execute("UPDATE users SET
                              first_name=split_part(name, ' ', 1)
                              last_name=split_part(name, ' ', 2);")

    alter_table(:users) do
      drop_column :name
    end
  end
end
