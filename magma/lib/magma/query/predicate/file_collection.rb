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
        require 'pry'
        binding.pry
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

        table.first[column_name] ? table.first[column_name].map do |f|
          binding.pry
          Magma.instance.storage.download_url(
            @model.project_name,
            f["filename"]
          ).to_s
        end : nil
      end
    end

    verb '::path' do
      child String

      extract do |table, identity|
        table.first[column_name] ? table.first[column_name].map { |f| f["filename"] } : nil
      end
    end

    verb '::original_filename' do
      child String

      extract do |table, identity|
        table.first[column_name] ? table.first[column_name].map { |f| f["original_filename"] } : nil
      end
    end

    verb '::md5' do
      child String

      extract do |table, identity|
        table.first[column_name]&.map { |f| @md5_set << f["filename"] }
      end
    end

    verb '::updated_at' do
      child String

      extract do |table, identity|
        table.first[column_name]&.map { |f| @updated_at_set << f["filename"] }
      end
    end

    verb '::all' do
      child String

      extract do |table, identity|
        table.first[column_name] ? table.first[column_name].map do |f|
          f.symbolize_keys.update({
            url: Magma.instance.storage.download_url(
              @model.project_name,
              f["filename"]
            ).to_s
          })
        end : nil
      end
    end

    def select
      @is_subselect ? [] : [ Sequel[alias_name][@column_name].as(column_name) ]
    end
  end
end
