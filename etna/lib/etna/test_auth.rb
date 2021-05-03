require_relative 'user'

# This is an authentication layer you can use in testing. It will make an
# Etna::User as usual that your controller can respond to; you can pass
# permissions for this user directly in the Authorization header
module Etna
  class TestAuth < Auth
    def self.token_header(params)
      token = Base64.strict_encode64(params.to_json)
      return [ 'Authorization', "Etna something.#{token}" ]
    end

    def self.token_param(params)
      token = Base64.strict_encode64(params.to_json)
      return [ Etna::Auth.etna_url_param(:authorization).to_s, "Etna #{token}" ]
    end

    def self.hmac_header(signature)
      return [ Etna::Auth.etna_url_param(:signature).to_s, signature ]
    end

    def self.hmac_params(params)
      return {
        expiration: params.delete(:expiration) || DateTime.now.iso8601,
        id: params.delete(:id) || 'etna',
        nonce: 'nonce',
        signature: params.delete(:signature) || 'invalid',
        headers: params.keys.join(',')
      }.merge(params).map do |item, value|
        [ Etna::Auth.etna_url_param(item).to_s, value ]
      end.to_h
    end

    def approve_user(request)
      token = auth(request,:etna)

      return false unless token

      # Here we simply base64-encode our user hash and pass it through
      # In order to behave more like "real" tokens, we expect the user hash to be
      #   in index 1 after splitting by ".".
      # We do this to support Metis client tests, we pass in tokens with multiple "."-separated parts, so
      #   have to account for that.
      payload = JSON.parse(Base64.decode64(token.split('.')[1]))
      request.env['etna.user'] = Etna::User.new(payload.map{|k,v| [k.to_sym, v]}.to_h, token)
    end

    def approve_hmac(request)
      hmac_signature = etna_param(request, :signature)

      return false unless hmac_signature

      headers = (etna_param(request, :headers)&.split(/,/) || []).map do |header|
        [ header.to_sym, etna_param(request, header) ]
      end.to_h

      hmac_params = {
        method: request.request_method,
        host: request.host,
        path: request.path,

        expiration: etna_param(request, :expiration) || DateTime.now.iso8601,
        id: etna_param(request, :id) || 'etna',
        nonce: etna_param(request, :nonce) || 'nonce',
        headers: headers,
        test_signature: hmac_signature
      }

      hmac = Etna::TestHmac.new(application, hmac_params)

      request.env['etna.hmac'] = hmac

      return nil unless hmac.valid?

      params(request).update(headers)

      return true
    end
  end
  class TestHmac < Hmac
    def valid?
      @test_signature == 'valid' || super
    end
  end
end
