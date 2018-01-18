module Etna
  class Hmac
    def initialize signer, params
      @signer = signer

      # These are the items that need to be signed
      [ :method, :host, :path, :timestamp, :nonce, :id, :headers ].each do |item|
        raise ArgumentError, "Hmac requires param #{item}" unless params[item]
        instance_variable_set("@#{item}", params[item])
      end

      raise ArgumentError, "Headers must be a Hash" unless @headers.is_a?(Hash)
    end

    def has_header? header, value
      @headers[header] == value
    end

    def valid_signature? test_signature
      signature == test_signature
    end

    def signature
      @signer.hmac(text_to_sign)
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
        @timestamp,
      ].join("\n")
    end
  end
end
