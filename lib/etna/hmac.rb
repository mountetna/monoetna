module Etna
  class Hmac
    def initialize application, params
      @application = application

      # These are the items that need to be signed
      [ :method, :host, :path, :expiration, :nonce, :id, :headers ].each do |item|
        raise ArgumentError, "Hmac requires param #{item}" unless params[item]
        instance_variable_set("@#{item}", params[item])
      end

      @id = @id.to_sym

      raise ArgumentError, "Invalid id #{@id}" unless @application.config(:hmac_keys).key?(@id)

      raise ArgumentError, "Headers must be a Hash" unless @headers.is_a?(Hash)
    end

    # this returns arguments for URI::HTTP.build
    def url_params
      params = {
        authorization: "Hmac #{signature}",
        expiration: @expiration,
        nonce: @nonce,
        id: @id.to_s,
        headers: @headers.keys.join(','),
      }.merge(@headers).map do |name, value|
        [
          "X-Etna-#{ name.to_s.split(/_/).map(&:capitalize).join('-') }",
          value
        ]
      end.to_h

      return {
        host: @host,
        path: @path,
        query: URI.encode_www_form(params)
      }
    end

    def valid_signature? test_signature
      signature == test_signature && DateTime.parse(@expiration) >= DateTime.now
    end

    def signature
      @application.sign.hmac(text_to_sign, @application.config(:hmac_keys)[@id])
    end

    private

    # This scheme is adapted from the Hawk spec
    # (github:hueniverse/hawk) and the Acquia HTTP
    # Hmac spec (github:acquia/http-hmac-spec)


    def text_to_sign
      [
        # these come from the route
        @method,
        @host,
        @path,

        # these are set as headers or params
        @nonce,
        @id,
        @headers.map{|l| l.join('=')}.join(';'),
        @expiration,
      ].join("\n")
    end
  end
end
