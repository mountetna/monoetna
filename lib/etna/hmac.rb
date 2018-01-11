def get_headers(param, pattern)
  param.select do |key,value|
    key =~ /^#{pattern}[\w-]+$/
  end.map do |key, value|
    [
      key.to_s.sub(/^#{pattern}/,'').downcase.to_sym, value
    ]
  end
end

module Etna
  class Hmac
    class AuthHeader < Hash
      SEP1 = ';'
      SEP2 = '='
      SEP3 = ','

      REQUIRED_AUTH_ITEMS = [
        :nonce,
        :id,
        :headers,
        :signature
      ]

      def initialize headers
        if headers.is_a?(Hash)
          @auth = headers
        elsif headers.is_a?(String)
          @auth = Hash[
            headers.split(/#{SEP1}/).map do |field|
              k, v = field.split(/#{SEP2}/)
              [ k.to_sym, v ]
            end
          ]
        else
          @auth = {}
        end
      end

      def key? key
        @auth.key? key
      end

      def [] key
        @auth[key]
      end

      def to_s
        @auth.map do |k,v|
          [ k, v ].join(SEP2)
        end.join(SEP1)
      end

      def required_headers
        @required_headers ||= @auth[:headers] ? @auth[:headers].split(/#{SEP3}/).map(&:to_sym) : []
      end

      def valid?
        REQUIRED_AUTH_ITEMS.all? {|i| @auth.key?(i) }
      end
    end

    def self.from_request(request, signer)
      auth = Etna::Hmac::AuthHeader.new(
        request.env['HTTP_AUTHORIZATION'][/\AHmac (.*)\z/,1]
      )

      # the authorization header must be in order
      return [] unless auth.valid?
      
      timestamp = request.env['HTTP_X_AUTHORIZATION_TIMESTAMP']

      # we must have the authorization timestamp
      return [] unless timestamp

      etna_headers = Hash[
        # If there is an HTTP header set, use this
        get_headers(request.env, 'HTTP_X_ETNA_') +
        # If not, use a parameter hash
        get_headers(request.env['rack.request.params'], 'X-Etna-')
      ]

      headers = auth.required_headers.map do |header|
        [ header, etna_headers[header] ]
      end.to_h

      # we must have all of the required headers (or parameters)
      return [] unless headers.all? {|h,v| !v.nil? }

      hmac = new(
        signer,
        method: request.request_method,
        host: request.host,
        path: request.path,
        nonce: auth[:nonce],
        headers: headers,
        timestamp: timestamp,
        id: auth[:id]
      )

      return [ hmac, auth[:signature] ]
    end

    def initialize signer, params
      @signer = signer

      # These are the items that need to be signed
      
      # items that come from the request
      @method = params[:method]
      @host = params[:host]
      @path = params[:path]
      @timestamp = params[:timestamp]

      # items in the auth header
      @nonce = params[:nonce]
      @id = params[:id]

      # headers is a Hash of parameter and values, 
      # taken from the headers or URL parameters
      @headers = params[:headers]
    end

    def valid_signature? sig
      signature = @signer.hmac(text_to_sign)

      return signature == @auth[:signature]
    end

    private

    # This scheme is adapted from the Hawk spec
    # (github:hueniverse/hawk) and the Acquia HTTP
    # Hmac spec (github:acquia/http-hmac-spec)


    def text_to_sign
      [
        @method,
        @host,
        @path,
        @nonce,
        @id,
        @headers.map{|l| l.join('=')}.join(';'),
        @timestamp,
      ].join("\n")
    end

    def authorization_timestamp
    end


    
  end
end
