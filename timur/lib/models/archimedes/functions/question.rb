module Archimedes
  class Question < Archimedes::Function
    def call
      query, *_ = @args

      host = Timur.instance.config(:magma).fetch(:host)

      client = Etna::Client.new(
        host,
        @token)
      
      query_params = {
        project_name: @project_name,
        query: query.to_values
      }

      query_answer = client.query(query_params)

      # Loop the data and set the data types returned from Magma.
      recursive_parse(query_answer[:answer]) do |item|
        case item
        when /^(\d{4})-(\d{2})-(\d{2})T(\d{2})\:(\d{2})\:(\d{2})[+-](\d{2})\:(\d{2})/
          DateTime.parse(item)
        when Array
          item.all?{|v| v.is_a?(Array) && v.length == 2}? Vector.new(item) : item
        else
          item
        end
      end
    end
    
    def recursive_parse(result, &block)
      case result
      when Array
        block.call(result.map{ |item| recursive_parse(item, &block) })
      when Hash
        block.call(Hash[result.map{ |key,item| [ key, recursive_parse(item, &block)] }])
      else
        block.call(result)
      end
    end
  end
end
