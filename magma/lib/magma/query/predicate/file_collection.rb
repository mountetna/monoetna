require 'json'

class Magma
  class FileCollectionPredicate < Magma::ColumnPredicate
    def initialize question, model, alias_name, attribute, *query_args
      super
      @md5_set = Md5Set.new(@question.user, @model)
      @updated_at_set = UpdatedAtSet.new(@question.user, @model)
    end

    verb '::url' do
      child String

      extract do |table, identity|
        # Problem is that ... table.first[column_name] can be nested
        #   So we have to dig into the lowest level where actual
        #   file_collection data is located, and then grab the
        #   download URLs for those.
        # And then make sure to keep the nesting...!
        #         [7] pry(#<Magma::FileCollectionPredicate>)> table
        # => [{:mmjdrnaehy_name=>"The Twelve Labors of Hercules",
        #   :jwitubrgwz_certificates=>
        #    [["Nemean Lion",
        #      [["Nemean Lion",
        #        [{"filename"=>"monster-Nemean Lion-certificates-0.txt", "original_filename"=>"sb_diploma_lion.txt"},
        #         {"filename"=>"monster-Nemean Lion-certificates-1.txt", "original_filename"=>"sm_diploma_lion.txt"}]]]],
        #     ["Lernean Hydra",
        #      [["Lernean Hydra",
        #        [{"filename"=>"monster-Lernean Hydra-certificates-0.txt", "original_filename"=>"ba_diploma_hydra.txt"},
        #         {"filename"=>"monster-Lernean Hydra-certificates-1.txt", "original_filename"=>"phd_diploma_hydra.txt"}]]]],
        #     ["Ceryneian Hind", [["Ceryneian Hind", nil]]]]}]
        if !table.first[column_name]
          nil
        elsif @is_subselect
          nested_reduce_and_apply(
            table.first[column_name],
            0,
            &request_download_url
          )
        else
          table.first[column_name].map do |f|
            Magma.instance.storage.download_url(
              @model.project_name,
              f["filename"]
            ).to_s
          end
        end
      end
    end

    verb '::path' do
      child String

      extract do |table, identity|
        if !table.first[column_name]
          nil
        elsif @is_subselect
          nested_reduce_and_apply(
            table.first[column_name],
            0,
            &extract_attribute("filename")
          )
        else
          table.first[column_name].map { |f| f["filename"] }
        end
      end
    end

    verb '::original_filename' do
      child String

      extract do |table, identity|
        if !table.first[column_name]
          nil
        elsif @is_subselect
          nested_reduce_and_apply(
            table.first[column_name],
            0,
            &extract_attribute("original_filename")
          )
        else
          table.first[column_name].map { |f| f["original_filename"] }
        end
      end
    end

    verb '::md5' do
      child String

      extract do |table, identity|
        if !table.first[column_name]
          nil
        elsif @is_subselect
          nested_reduce_and_apply(
            table.first[column_name],
            0,
            &append_attribute_to_set(@md5_set, "filename")
          )
        else
          table.first[column_name].map { |f| @md5_set << f["filename"] }
        end
      end
    end

    verb '::updated_at' do
      child String

      extract do |table, identity|
        if !table.first[column_name]
          nil
        elsif @is_subselect
          nested_reduce_and_apply(
            table.first[column_name],
            0,
            &append_attribute_to_set(@updated_at_set, "filename")
          )
        else
          table.first[column_name].map { |f| @updated_at_set << f["filename"] }
        end
      end
    end

    verb '::all' do
      child String

      extract do |table, identity|
        if !table.first[column_name]
          nil
        elsif @is_subselect
          nested_reduce_and_apply(
            table.first[column_name],
            0,
            &append_attribute_to_set(@updated_at_set, "filename")
          )
        else
          table.first[column_name].map do |f|
            f.symbolize_keys.update({
              url: request_download_url.call(f)
            })
          end
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

    def nested_reduce_and_apply(nested_array_data, level, &file_transformation)
      nested_array_data&.reduce([]) do |result, value|
        if include_identifier?(value, level)
          # is an tuple, [identifier, data]
          identifier = value.first
          data = nested_reduce_and_apply(value.last, level + 1, &file_transformation).flatten

          result << [identifier, data]
        elsif skip_identifier?(value, level)
          # Skip the intervening identifiers
          result << nested_reduce_and_apply(value.last, level + 1, &file_transformation)
        else
          result << file_transformation.call(value)
        end
      end
    end

    def include_identifier?(value, level)
      value.is_a?(Array) && 0 == level
    end

    def skip_identifier?(value, level)
      value.is_a?(Array) && level > 0
    end
  end
end
