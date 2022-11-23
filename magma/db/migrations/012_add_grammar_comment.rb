Sequel.migration do
  up do
    alter_table(:grammars) do
      add_column :comment, String
    end
    
    Magma.instance.db.execute("UPDATE grammars SET comment='update' WHERE grammars.comment is NULL")

    alter_table(:grammars) do
      set_column_not_null :comment
    end
  end

  down do
    alter_table(:grammars) do
      drop_column :comment
    end
  end
end
