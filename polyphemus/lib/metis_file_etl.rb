require 'ostruct'

class Polyphemus
  # Abstract base class meant to do things.
  class MetisFileEtl
    attr_reader :metis_client, :cursor_group

    def initialize(metis_client:, cursor_group:)
      @metis_client = metis_client
      @cursor_group = cursor_group
    end

    # Subclasses should override if they wish to adjust or add to the params of the find request.
    def prepare_find_request(find_request)
    end

    def execute
      cursor_group.with_next do |cursor|
        find_request = Etna::Clients::Metis::FindRequest.new(
            project_name: cursor[:project_name],
            bucket_name: cursor[:bucket_name],
            limit: cursor[:limit] || 20,
            offset: cursor[:offset] || 0,
            params: [
                Etna::Clients::Metis::FindParam.new(
                    type: 'file',
                    attribute: 'updated_at',
                    predicate: '>=',
                    value: cursor[:updated_at],
                )
            ]
        )

        prepare_find_request(find_request)
        response = metis_client.find(find_request)
        process(cursor, response.files)
        update_cursor(cursor, find_request, response)
      end
    end

    # Subclasses should override with their processing loop.
    # IMPORTANT!  Keep in mind that because the process can fail part way through work, there is no guarantee of single
    # delivery and that process should be idempotent to the best of its ability.
    # This class does best effort to manage work 'split' along the same updated at to try to ensure that for any
    # successful process and update_cursor call, a duplicate one should not occur after.
    def process(cursor, files)
    end

    def update_cursor(cursor, find_request, response)
      new_updated_at = response.files.all.map(&:updated_at).max
      ids = response.files.all.map(&:id)

      if new_updated_at.nil?
        return
      end

      if new_updated_at == cursor[:updated_at]
        cursor[:ids] += ids
      else
        cursor[:ids] = ids
      end

      cursor[:updated_at] = new_updated_at
    end
  end
end
