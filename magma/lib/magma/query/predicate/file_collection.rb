require 'json'
require_relative 'with_aggregated_data_module'

class Magma
  class FileCollectionPredicate < Magma::ColumnPredicate
    include WithAggregatedDataModule

    def initialize question, model, alias_name, attribute, *query_args
      super
      @md5_set = Md5Set.new(@question.user, @model)
      @updated_at_set = UpdatedAtSet.new(@question.user, @model)
    end

    verb '::url' do
      child String

      extract do |table, identity|
        if !table.first[column_name]
          Magma::NilAnswer.new
        elsif @is_subselect
          nested_reduce_and_apply(
            table.first[column_name],
            0,
            &request_download_url
          )
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
        elsif @is_subselect
          nested_reduce_and_apply(
            table.first[column_name],
            0,
            &extract_attribute("filename")
          )
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
        elsif @is_subselect
          nested_reduce_and_apply(
            table.first[column_name],
            0,
            &extract_attribute("original_filename")
          )
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
        elsif @is_subselect
          nested_reduce_and_apply(
            table.first[column_name],
            0,
            &append_attribute_to_set(@md5_set, "filename")
          )
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
        elsif @is_subselect
          nested_reduce_and_apply(
            table.first[column_name],
            0,
            &append_attribute_to_set(@updated_at_set, "filename")
          )
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
        elsif @is_subselect
          nested_reduce_and_apply(
            table.first[column_name],
            0,
            &append_attribute_to_set(@updated_at_set, "filename")
          )
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

    def select
      @is_subselect ? [] : [ Sequel[alias_name][@column_name].as(column_name) ]
    end

    private

    def update_download_url
      lambda { |file_hash|
        file_hash.symbolize_keys.update({
          url: request_download_url.call(file_hash)
        })
      }
    end

    def request_download_url
      lambda { |file_hash|
        Magma.instance.storage.download_url(
          @model.project_name,
          file_hash["filename"]
        ).to_s
      }
    end

    def extract_attribute(file_attribute_name)
      lambda { |file_hash|
        file_hash[file_attribute_name]
      }
    end

    def append_attribute_to_set(metis_set, file_attribute_name)
      lambda { |file_hash|
        metis_set << file_hash[file_attribute_name]
      }
    end

    def as_answer_array(raw_answers)
      raw_answers.map do |answer|
        Magma::Answer.new(answer)
      end
    end
  end
end
