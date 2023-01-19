require_relative 'column_subselect'
require_relative '../with_file_collection_predicate_methods'

class Magma
  class FileCollectionSubselectPredicate < Magma::ColumnSubselectPredicate
    include WithFileCollectionPredicateMethods

    def self.verbs
      Magma::FileCollectionPredicate.verbs.merge(@verbs)
    end

    verb '::url' do
      child String

      extract do |table, identity|
        if !table.first[column_name]
          Magma::NilAnswer.new
        else
          nested_reduce_and_apply(
            table.first[column_name],
            0,
            &request_download_url
          )
        end
      end
    end

    verb '::path' do
      child String

      extract do |table, identity|
        if !table.first[column_name]
          Magma::NilAnswer.new
        else
          nested_reduce_and_apply(
            table.first[column_name],
            0,
            &extract_attribute("filename")
          )
        end
      end
    end

    verb '::original_filename' do
      child String

      extract do |table, identity|
        if !table.first[column_name]
          Magma::NilAnswer.new
        else
          nested_reduce_and_apply(
            table.first[column_name],
            0,
            &extract_attribute("original_filename")
          )
        end
      end
    end

    verb '::md5' do
      child String

      extract do |table, identity|
        if !table.first[column_name]
          Magma::NilAnswer.new
        else
          nested_reduce_and_apply(
            table.first[column_name],
            0,
            &append_attribute_to_set(@md5_set, "filename")
          )
        end
      end
    end

    verb '::updated_at' do
      child String

      extract do |table, identity|
        if !table.first[column_name]
          Magma::NilAnswer.new
        else
          nested_reduce_and_apply(
            table.first[column_name],
            0,
            &append_attribute_to_set(@updated_at_set, "filename")
          )
        end
      end
    end

    verb '::all' do
      child String

      extract do |table, identity|
        if !table.first[column_name]
          Magma::NilAnswer.new
        else
          nested_reduce_and_apply(
            table.first[column_name],
            0,
            &append_attribute_to_set(@updated_at_set, "filename")
          )
        end
      end
    end

    def select
      []
    end
  end
end
