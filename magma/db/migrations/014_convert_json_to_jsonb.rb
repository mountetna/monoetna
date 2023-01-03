Sequel.migration do
  class Helper

    def self.convert_attribute_columns(migration, new_type)
      matching_model_attributes.each do |table_name_tuple, attribute_names|
        migration.alter_table(Sequel[table_name_tuple.first][table_name_tuple.last]) do
          attribute_names.each do |attribute_name|
            set_column_type attribute_name, new_type, using: Sequel.cast(attribute_name, new_type)
          end
        end
      end
    end

    private

    def self.matching_model_attributes
      {}.tap do |table_attributes|
        Magma.instance.db[:attributes].where(type: json_type_attributes).each do |attribute|
          schema_name = attribute[:project_name].to_sym
          table_name = attribute[:model_name].pluralize.to_sym

          table_tuple = [schema_name, table_name]

          table_attributes[table_tuple] = [] unless table_attributes.keys.include?(table_tuple)

          table_attributes[table_tuple] << attribute[:attribute_name].to_sym
        end
      end
    end

    def self.json_type_attributes
      [ 'file', 'image', 'file_collection', 'matrix' ]
    end
  end

  up do
    new_type = :jsonb

    alter_table(:grammars) do
      set_column_type :config, new_type, using: Sequel.cast(:config, new_type)
    end
    alter_table(:attributes) do
      set_column_type :validation, new_type, using: Sequel.cast(:validation, new_type)
    end
    alter_table(:models) do
      set_column_type :dictionary, new_type, using: Sequel.cast(:dictionary, new_type)
    end

    Helper.convert_attribute_columns(self, new_type)
  end

  down do
    new_type = :json

    alter_table(:grammars) do
      set_column_type :config, new_type, using: Sequel.cast(:config, new_type)
    end
    alter_table(:attributes) do
      set_column_type :validation, new_type, using: Sequel.cast(:validation, new_type)
    end
    alter_table(:models) do
      set_column_type :dictionary, new_type, using: Sequel.cast(:dictionary, new_type)
    end

    Helper.convert_attribute_columns(self, new_type)
  end
end
