require 'webmock/rspec'
require 'vcr'
require 'openssl'
require 'digest/sha2'
require 'base64'

def setup_base_vcr(spec_helper_dir)
  VCR.configure do |c|
    c.hook_into :webmock
    c.cassette_serializers
    c.cassette_library_dir = ::File.join(spec_helper_dir, 'fixtures', 'cassettes')
    c.allow_http_connections_when_no_cassette = true

    c.register_request_matcher :try_body do |request_1, request_2|
      if request_1.headers['Content-Type'].first =~ /application\/json/
        if request_2.headers['Content-Type'].first =~ /application\/json/
          request_1_json = begin JSON.parse(request_1.body) rescue 'not-json' end
          request_2_json = begin JSON.parse(request_2.body) rescue 'not-json' end
          request_1_json == request_2_json
        else
          false
        end
      else
        request_1.body == request_2.body
      end
    end

    c.default_cassette_options = {
        serialize_with: :compressed,
        record: if ENV['IS_CI'] == '1'
          :none
        else
          ENV['RERECORD'] ? :all : :once
        end,
        match_requests_on: [:method, :uri, :try_body]
    }

    # Filter the authorization headers of any request by replacing any occurrence of that request's
    # Authorization value with <AUTHORIZATION>
    c.filter_sensitive_data('<AUTHORIZATION>') do |interaction|
      interaction.request.headers['Authorization'].first
    end

    c.before_record do |interaction|
      key = prepare_vcr_secret

      cipher = OpenSSL::Cipher.new("AES-256-CBC")
      iv = cipher.random_iv

      cipher.encrypt
      cipher.key = key
      cipher.iv = iv

      encrypted = cipher.update(interaction.response.body)
      encrypted << cipher.final

      interaction.response.body = [iv, encrypted].pack('mm')
    end

    c.before_playback do |interaction|
      key = prepare_vcr_secret

      iv, encrypted = interaction.response.body.unpack('mm')

      cipher = OpenSSL::Cipher.new("AES-256-CBC")
      cipher.decrypt
      cipher.key = key
      cipher.iv = iv

      plain = cipher.update(encrypted)
      plain << cipher.final

      interaction.response.body = plain
    end
  end
end

def prepare_vcr_secret
  secret = ENV["CI_SECRET"]

  if (secret.nil? || secret.empty?) && ENV['IS_CI'] != '1'
    current_example = RSpec.current_example
    RSpec::Core::Pending.mark_pending! current_example, 'CI_SECRET must be set to run this test'
    raise "CI_SECRET must be set to run this test"
  end

  digest = Digest::SHA256.new
  digest.update(secret)
  digest.digest
end