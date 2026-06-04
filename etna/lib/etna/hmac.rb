module Etna
  class Hmac
    # These are the items that need to be signed
    SIGN_ITEMS=[ :method, :host, :path, :expiration, :nonce, :id, :headers ]
    SIGN_ITEMS.each { |item| attr_reader item }

    def initialize application, params
      @application = application

      @test_signature = params.delete(:test_signature)

      SIGN_ITEMS.each do |item|
        raise ArgumentError, "Hmac requires param #{item}" unless params[item]
        instance_variable_set("@#{item}", params[item])
      end

      @id = @id.to_sym

      raise ArgumentError, "Headers must be a Hash" unless @headers.is_a?(Hash)
    end

    # this returns arguments for URI::HTTP.build
    def url_params(with_headers=true)
      params = {
        signature: signature,
        expiration: @expiration,
        nonce: @nonce,
        id: @id.to_s,
        headers: @headers.keys.join(','),
      }.merge(with_headers ? @headers : {}).map do |name, value|
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

    def valid?
      valid_id? && valid_signature? && valid_timestamp?
    end

    def signature
      @application.sign.hmac(text_to_sign, @application.config(@id)[:hmac_key])
    end

    private

    def valid_signature?
      signature == @test_signature
    end

    def valid_timestamp?
      DateTime.parse(@expiration) >= DateTime.now
    end

    def valid_id?
      return false if !@application.config(@id)
      @application.config(@id).key?(:hmac_key)
    end

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
        JSON.generate(@headers),
        @expiration,
      ].join("\n")
    end
  end
end
