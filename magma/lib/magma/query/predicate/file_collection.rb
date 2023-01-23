require 'json'
require_relative 'with_file_collection_predicate_methods'

class Magma
  class FileCollectionPredicate < Magma::ColumnPredicate
    include WithFileCollectionPredicateMethods

    verb '::url' do
      child String

      extract do |table, identity|
        if !table.first[column_name]
          Magma::NilAnswer.new
        else
          as_answer_array(
            table.first[column_name].map do |f|
              request_download_url.call(f)
            end
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
          as_answer_array(
            table.first[column_name].map { |f| f["filename"] }
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
          as_answer_array(
            table.first[column_name].map { |f| f["original_filename"] }
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
          as_answer_array(
            table.first[column_name].map { |f| @md5_set << f["filename"] }
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
          as_answer_array(
            table.first[column_name].map { |f| @updated_at_set << f["filename"] }
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
          as_answer_array(
            table.first[column_name].map do |f|
              f.symbolize_keys.update({
                url: request_download_url.call(f)
              })
            end
          )
        end
      end
    end

    def select(incoming_alias_name=nil, incoming_attribute=nil)
      [ Sequel[alias_name][@column_name].as(column_name) ]
    end

    private

    def as_answer_array(raw_answers)
      raw_answers.map do |answer|
        Magma::Answer.new(answer)
      end
    end
  end
end
