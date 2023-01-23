class Magma
  module WithFileCollectionPredicateMethods
    def initialize question, model, alias_name, attribute, *query_args
      super
      @md5_set = Md5Set.new(@question.user, @model)
      @updated_at_set = UpdatedAtSet.new(@question.user, @model)
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
  end
end