# General signing/hashing utilities.
require 'jwt'
require 'securerandom'

module Etna
  class SignService
    def initialize(application)
      @application = application
    end

    def hash_password(password)
      signature(
        [ password, @application.config(:pass_salt) ],
        @application.config(:pass_algo)
      )
    end

    def hmac(message, key)
      OpenSSL::HMAC.hexdigest(
        'SHA256',
        key,
        message
      )
    end

    def jwt_token(payload)
      return JWT.encode(
        payload,
        private_key,
        @application.config(:token_algo)
      )
    end

    def uid(size=nil)
      SecureRandom.hex(size)
    end

    def jwt_decode(token)
      return JWT.decode(
        token,
        public_key,
        true,
        algorithm: @application.config(:token_algo)
      )
    end

    def private_key
      @private_key ||= OpenSSL::PKey::RSA.new(@application.config(:rsa_private))
    end

    def public_key
      @public_key ||= OpenSSL::PKey::RSA.new(@application.config(:rsa_public))
    end

    def generate_private_key(key_size)
      OpenSSL::PKey::RSA.generate(key_size)
    end

    private

    def signature(params, algo)
      algo = algo.upcase.to_sym
      raise "Unknown signature algorithm!" unless [ :MD5, :SHA256 ].include?(algo)
      Digest.const_get(algo).hexdigest(params.join)
    end
  end
end
