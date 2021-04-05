require 'openssl'

class Janus
  class Nonce
    def self.nonce_key
      @nonce_key ||= OpenSSL::PKey::RSA.new(2048)
    end

    def self.valid_nonce?(nonce)
      timestamp, nonce_sig = Base64.decode64(nonce).split(/\./)

      return false unless timestamp && nonce_sig

      begin
        date = DateTime.parse(timestamp)
      rescue ArgumentError
        # Invalid date
        return false
      end

      return false if (DateTime.now - date) * 24 * 60 * 60 > 60

      return Janus::Nonce.new(timestamp).to_s == nonce
    end

    def initialize timestamp
      @timestamp = timestamp
    end

    def to_s
      # sign the time with our nonce key
      signature = Janus::Nonce.nonce_key.sign(
        OpenSSL::Digest::SHA256.new,
        @timestamp
      )

      # Digest the signature
      hash = Digest::SHA256.hexdigest(signature)

      # encode the pair
      return Base64.strict_encode64("#{@timestamp}.#{hash}")
    end
  end
end
