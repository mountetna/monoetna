describe Magma::Migration do
  before(:each) { Timecop.freeze('2000-01-01') }
  after(:each) { Timecop.return }

  context 'empty migrations' do
    it 'does nothing if there is no change' do
      migration = Labors::Project.migration

      expect(migration).to be_empty
      expect(migration.to_s).to eq('')
    end
  end

  context 'creation migrations' do
    after(:each) do
      Labors.send(:remove_const, :Olympian)
      project = Magma.instance.get_project(:labors)
      project.models.delete(:olympian)
    end

    it 'suggests a creation migration for identifiers' do
      module Labors
        class Olympian < Magma::Model
          identifier :name, type: String
        end
      end
      migration = Labors::Olympian.migration
      expect(migration.to_s).to eq <<EOT.chomp
    create_table(Sequel[:labors][:olympians]) do
      primary_key :id
      DateTime :created_at
      DateTime :updated_at
      String :name
      unique :name
    end
EOT
    end

    it 'suggests a creation migration for attributes' do
      module Labors
        class Olympian < Magma::Model
          integer :number
        end
      end
      migration = Labors::Olympian.migration
      expect(migration.to_s).to eq <<EOT.chomp
    create_table(Sequel[:labors][:olympians]) do
      primary_key :id
      DateTime :created_at
      DateTime :updated_at
      Integer :number
    end
EOT
    end

    it 'suggests nothing in the creation for a table or collection attribute' do
      module Labors
        class Olympian < Magma::Model
          collection :victim
          table :prize
        end
      end
      migration = Labors::Olympian.migration
      expect(migration.to_s).to eq <<EOT.chomp
    create_table(Sequel[:labors][:olympians]) do
      primary_key :id
      DateTime :created_at
      DateTime :updated_at
    end
EOT
    end

    it 'suggests a creation migration for json attributes' do
      module Labors
        class Olympian < Magma::Model
          match :prayers
        end
      end
      migration = Labors::Olympian.migration
      expect(migration.to_s).to eq <<EOT.chomp
    create_table(Sequel[:labors][:olympians]) do
      primary_key :id
      DateTime :created_at
      DateTime :updated_at
      json :prayers
    end
EOT
    end

    it 'suggests a creation migration for link attributes' do
      module Labors
        class Olympian < Magma::Model
          parent :project

          link :monster
        end
      end
      migration = Labors::Olympian.migration
      expect(migration.to_s).to eq <<EOT.chomp
    create_table(Sequel[:labors][:olympians]) do
      primary_key :id
      DateTime :created_at
      DateTime :updated_at
      foreign_key :project_id, Sequel[:labors][:projects]
      index :project_id
      foreign_key :monster_id, Sequel[:labors][:monsters]
      index :monster_id
    end
EOT
    end
  end

  context 'update migrations' do
    def remove_attribute(model, attribute)
      model.attributes.delete(attribute)
      model.instance_variable_set("@identity",nil) if model.identity.column_name.to_sym == attribute
    end

    def remove_link(model, attribute_name)
      attribute = model.attributes[attribute_name]
      reciprocal_model = Magma.instance.get_model('labors', attribute.link_model_name)
      reciprocal_model.attributes.delete(attribute.link_attribute_name.to_sym)
      model.attributes.delete(attribute_name)
      model.instance_variable_set("@identity",nil) if model.identity.column_name.to_sym == attribute_name
    end

    it 'suggests an update migration for identifiers' do
      module Labors
        class Prize < Magma::Model
          identifier :prize_code, type: String
        end
      end
      migration = Labors::Prize.migration
      expect(migration.to_s).to eq <<EOT.chomp
    alter_table(Sequel[:labors][:prizes]) do
      add_column :prize_code, String
      add_unique_constraint :prize_code
    end
EOT
      remove_attribute(Labors::Prize, :prize_code)
      Labors::Prize.order()
    end

    it 'suggests an update migration for attributes' do
      module Labors
        class Prize < Magma::Model
          float :weight
        end
      end
      migration = Labors::Prize.migration
      expect(migration.to_s).to eq <<EOT.chomp
    alter_table(Sequel[:labors][:prizes]) do
      add_column :weight, Float
    end
EOT
      remove_attribute(Labors::Prize, :weight)
    end

    it 'suggests an update migration for json attributes' do
      module Labors
        class Prize < Magma::Model
          match :dimensions
        end
      end
      migration = Labors::Prize.migration
      expect(migration.to_s).to eq <<EOT.chomp
    alter_table(Sequel[:labors][:prizes]) do
      add_column :dimensions, :json
    end
EOT
      remove_attribute(Labors::Prize, :dimensions)
    end

    it 'suggests an update migration for link attributes' do
      module Labors
        class Prize < Magma::Model
          link :monster
        end
      end
      migration = Labors::Prize.migration
      expect(migration.to_s).to eq <<EOT.chomp
    alter_table(Sequel[:labors][:prizes]) do
      add_foreign_key :monster_id, Sequel[:labors][:monsters]
      add_index :monster_id
    end
EOT
      remove_attribute(Labors::Prize, :monster)
    end

    it 'removes attributes' do
      worth = Labors::Prize.attributes[:worth].dup
      remove_attribute(Labors::Prize,:worth)

      migration = Labors::Prize.migration
      expect(migration.to_s).to eq <<EOT.chomp
    alter_table(Sequel[:labors][:prizes]) do
      drop_foreign_constraint_if_exists :#{worth.column_name}
      rename_column :#{worth.column_name}, :#{worth.column_name}_946684800_backup
    end
EOT
      Labors::Prize.attributes[:worth] = worth
    end

    it 'removes attributes when column name does not match attribute name' do
      worth = Labors::Prize.attributes[:worth].dup

      action = Magma::RenameAttributeAction.new("labors", {
        action: "rename_attribute",
        model_name: "prize",
        attribute_name: "worth",
        new_attribute_name: "worthiness"
      })
      action.perform

      expect(Labors::Prize.attributes.keys.include?(:worthiness)).to eq(true)
      remove_attribute(Labors::Prize,:worthiness)

      migration = Labors::Prize.migration
      expect(migration.to_s).to eq <<EOT.chomp
    alter_table(Sequel[:labors][:prizes]) do
      drop_foreign_constraint_if_exists :#{worth.column_name}
      rename_column :#{worth.column_name}, :#{worth.column_name}_946684800_backup
    end
EOT
      Labors::Prize.attributes[:worth] = worth
    end

    xit 'makes no changes when removing child, collection or table attributes' do
      monster = Labors::Labor.attributes[:monster]
      prize = Labors::Labor.attributes[:prize]
      victim = Labors::Monster.attributes[:victim]
      remove_attribute(Labors::Labor,:monster)
      remove_attribute(Labors::Labor,:prize)
      remove_attribute(Labors::Monster,:victim)

      expect(Labors::Labor.migration).to be_empty
      expect(Labors::Monster.migration).to be_empty

      Labors::Labor.attributes[:monster] = monster
      Labors::Labor.attributes[:prize] = prize
      Labors::Monster.attributes[:victim] = victim
    end

    it 'removes link and reciprocal' do
      habitat = Labors::Monster.attributes[:habitat].dup
      monster = Labors::Habitat.attributes[:monster].dup
      remove_link(Labors::Monster,:habitat)

      migration = Labors::Monster.migration
      expect(migration.to_s).to eq <<EOT.chomp
    alter_table(Sequel[:labors][:monsters]) do
      drop_foreign_constraint_if_exists :#{habitat.column_name}
      rename_column :#{habitat.column_name}, :#{habitat.column_name}_946684800_backup
    end
EOT

      reciprocal_migration = Labors::Habitat.migration
      # Because the FK is on Labors::Monster, no migration needed
      #   for Labors::Habitat.
      expect(reciprocal_migration.to_s).to be_empty
      Labors::Monster.attributes[:habitat] = habitat
      Labors::Habitat.attributes[:monster] = monster
    end

    it 'removes link and reciprocal when calling remove_link on parent attribute' do
      habitat = Labors::Monster.attributes[:habitat].dup
      monster = Labors::Habitat.attributes[:monster].dup
      remove_link(Labors::Habitat,:monster)

      migration = Labors::Monster.migration
      expect(migration.to_s).to eq <<EOT.chomp
    alter_table(Sequel[:labors][:monsters]) do
      drop_foreign_constraint_if_exists :#{habitat.column_name}
      rename_column :#{habitat.column_name}, :#{habitat.column_name}_946684800_backup
    end
EOT

      reciprocal_migration = Labors::Habitat.migration
      # Because the FK is on Labors::Monster, no migration needed
      #   for Labors::Habitat.
      expect(reciprocal_migration.to_s).to be_empty
      Labors::Monster.attributes[:habitat] = habitat
      Labors::Habitat.attributes[:monster] = monster
    end

    it 'removes link and reciprocal when column name does not match attribute name' do
      habitat = Labors::Monster.attributes[:habitat].dup
      monster = Labors::Habitat.attributes[:monster].dup

      action = Magma::RenameAttributeAction.new("labors", {
        action: "rename_attribute",
        model_name: "monster",
        attribute_name: "habitat",
        new_attribute_name: "habitat_name"
      })
      action.perform

      expect(Labors::Monster.attributes.keys.include?(:habitat_name)).to eq(true)
      remove_link(Labors::Monster,:habitat_name)

      migration = Labors::Monster.migration
      expect(migration.to_s).to eq <<EOT.chomp
    alter_table(Sequel[:labors][:monsters]) do
      drop_foreign_constraint_if_exists :#{habitat.column_name}
      rename_column :#{habitat.column_name}, :#{habitat.column_name}_946684800_backup
    end
EOT

      reciprocal_migration = Labors::Habitat.migration
      # Because the FK is on Labors::Monster, no migration needed
      #   for Labors::Habitat.
      expect(reciprocal_migration.to_s).to be_empty
      Labors::Monster.attributes[:habitat] = habitat
      Labors::Habitat.attributes[:monster] = monster
    end

    it 'makes multiple changes at once' do
      worth = Labors::Prize.attributes[:worth]
      remove_attribute(Labors::Prize,:worth)
      module Labors
        class Prize < Magma::Model
          float :weight
        end
      end

      migration = Labors::Prize.migration
      expect(migration.to_s).to eq <<EOT.chomp
    alter_table(Sequel[:labors][:prizes]) do
      add_column :weight, Float
      drop_foreign_constraint_if_exists :#{worth.column_name}
      rename_column :#{worth.column_name}, :#{worth.column_name}_946684800_backup
    end
EOT
      remove_attribute(Labors::Prize, :weight)
      Labors::Prize.attributes[:worth] = worth
    end

    context 'changes column types' do
      it 'for simple types' do
        original_attribute = Labors::Prize.attributes.delete(:worth)
        Labors::Prize.attributes[:worth] = Magma::Model.float(
          :worth,
          column_name: original_attribute.column_name
        )

        migration = Labors::Prize.migration
        expect(migration.to_s).to eq <<EOT.chomp
    alter_table(Sequel[:labors][:prizes]) do
      set_column_type :#{original_attribute.column_name}, Float
    end
EOT

        Labors::Prize.attributes[:worth] = original_attribute
      end

      it 'for symbol types' do
        original_attribute = Labors::Prize.attributes.delete(:worth)
        Labors::Prize.attributes[:worth] = Magma::Model.image(
          :worth,
          column_name: original_attribute.column_name
        )

        migration = Labors::Prize.migration
        expect(migration.to_s).to eq <<EOT.chomp
    alter_table(Sequel[:labors][:prizes]) do
      set_column_type :#{original_attribute.column_name}, :json, using: 'worth::json'
    end
EOT
        Labors::Prize.attributes[:worth] = original_attribute
      end
    end
  end

  context 'remove model migrations' do
    let(:now) { DateTime.now.to_time.to_i }
    let(:backup_table_name) { "sidekicks_#{now}_backup" }

    before(:each) do
      Timecop.freeze('2000-01-01') # 946684800 since epoch
      @project = Magma.instance.get_project("labors")
      @sidekick_model = Magma.instance.db[:models].where(project_name: 'labors', model_name: 'sidekick').first
      @sidekick_attributes = @project.models[:sidekick].attributes.values
      @reciprocal_attribute = @project.models[:victim].attributes[:sidekick].dup
    end

    after(:each) do
      Timecop.return
      model = @project.load_model(@sidekick_model)
      model.load_attributes(@sidekick_attributes)
      @project.models[model.model_name] = model
      @project.models[:victim].load_attributes([@reciprocal_attribute])
    end

    it 'suggests a remove model migration and removes indices' do
      @project.models[:victim].attributes.delete(:sidekick)

      migration = @project.migrations.first

      expect(migration.to_s).to eq <<EOT.chomp
    alter_table(Sequel[:labors][:sidekicks]) do
      drop_index [:victim_id]
    end
    rename_table(Sequel[:labors][:sidekicks], Sequel[:labors][:#{backup_table_name}])
EOT
    end
  end
end
