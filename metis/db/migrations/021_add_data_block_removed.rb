Sequel.migration do
    change do
      alter_table(:data_blocks) do
        add_column :removed, TrueClass, default: false
      end
    end
  end
