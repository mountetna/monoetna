require_relative '../md5_set'
require_relative '../updated_at_set'

class Magma
  class FilePredicate < Magma::ColumnPredicate
    attr_reader :requested_file_paths

    def initialize question, model, alias_name, attribute, *query_args
      super
      @md5_set = MetisMetadata.new(@question.user, @model, :file_hash)
      @file_size_set = MetisMetadata.new(@question.user, @model, :size)
      @updated_at_set = MetisMetadata.new(@question.user, @model, :updated_at)
    end

    verb '::url' do
      child String

      extract do |table, identity|
        table.first[column_name] ?
          Magma.instance.storage.download_url(
            @model.project_name,
            table.first[column_name]["filename"]
          ) :
          Magma::NilAnswer.new
      end
    end


    class MetisMetadataValue
      def initialize(predicate, file)
        @predicate = predicate
        @file = file

        @predicate.requested_file_paths << file
      end
    end

    verb '::md5' do
      child String

      extract do |table, identity|
        table.first[column_name] ?
          @md5_set << table.first[column_name]["filename"] :
          Magma::NilAnswer.new
      end
    end

    verb '::size' do
      child String

      extract do |table, identity|
        table.first[column_name] ?
          @file_size_set << table.first[column_name]["filename"] :
          Magma::NilAnswer.new
      end
    end

    verb '::updated_at' do
      child String

      extract do |table, identity|
        table.first[column_name] ?
          @updated_at_set << table.first[column_name]["filename"] :
          Magma::NilAnswer.new
      end
    end

    verb '::path' do
      child String

      extract do |table, identity|
        table.first[column_name] ?
          Magma::Answer.new(table.first[column_name]["filename"]) :
          Magma::NilAnswer.new
      end
    end

    verb '::original_filename' do
      child String

      extract do |table, identity|
        table.first[column_name] ?
          Magma::Answer.new(table.first[column_name]["original_filename"]) :
          Magma::NilAnswer.new
      end
    end

    verb '::all' do
      child String

      extract do |table, identity|
        table.first[column_name] ?
          Magma::Answer.new(table.first[column_name].symbolize_keys) :
          Magma::NilAnswer.new
      end
    end

    verb [ '::equals' ], String do
      child TrueClass

      constraint do
        json_constraint(@column_name, "filename", @arguments[1])
      end
    end

    verb [ '::not' ], String do
      child TrueClass

      constraint do
        not_json_constraint(@column_name, "filename", @arguments[1])
      end
    end

    def select(incoming_alias_name=nil, incoming_attribute=nil)
      [ Sequel[alias_name][@column_name].as(column_name) ]
    end
  end
end
