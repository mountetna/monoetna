require 'cgi'
require 'webmock/rspec'
require 'vcr'
require 'openssl'
require 'digest/sha2'
require 'base64'

def clean_query(json_or_string)
  if json_or_string.is_a?(Hash) && json_or_string.include?('upload_path')
    json_or_string['upload_path'] = clean_query(json_or_string['upload_path'])
    json_or_string
  elsif json_or_string.is_a?(String)
    uri = URI(json_or_string)

    if uri.query&.include?('X-Etna-Signature')
      uri.query = 'etna-signature'
    end

    uri.to_s
  else
    json_or_string
  end
end

def setup_base_vcr(spec_helper_dir, server: nil, application: nil)
  VCR.configure do |c|
    c.hook_into :webmock
    c.cassette_serializers
    c.cassette_library_dir = ::File.join(spec_helper_dir, 'fixtures', 'cassettes')
    c.allow_http_connections_when_no_cassette = true

    c.register_request_matcher :verify_uri_route do |request_1, request_2|
      next true if server.nil? || application.nil?

      route_match = request_1.uri =~ /https:\/\/#{application.dev_route}(.*)/
      if route_match && route_match[1]
        def request_1.request_method
          method.to_s.upcase
        end

        def request_1.path
          URI(uri).path
        end

        !!server.find_route(request_1)
      else
        true
      end
    end

    c.register_request_matcher :try_uri do |request_1, request_2|
      clean_query(request_1.uri) == clean_query(request_2.uri)
    end

    c.register_request_matcher :try_body do |request_1, request_2|
      if request_1.headers['Content-Type'].first =~ /application\/json/
        if request_2.headers['Content-Type'].first =~ /application\/json/
          request_1_json = begin
            JSON.parse(request_1.body) rescue 'not-json'
          end
          request_2_json = begin
            JSON.parse(request_2.body) rescue 'not-json'
          end

          clean_query(request_1_json) == clean_query(request_2_json)
        else
          false
        end
      else
        request_1.body == request_2.body
      end
    end

    if File.exists?('log')
      c.debug_logger = File.open('log/vcr_debug.log', 'w')
    end

    c.default_cassette_options = {
        serialize_with: :compressed,
        record: if ENV['IS_CI'] == '1'
          :none
        else
          ENV['RERECORD'] ? :all : :once
        end,
        match_requests_on: [:method, :try_uri, :try_body, :verify_uri_route]
    }

    # Filter the authorization headers of any request by replacing any occurrence of that request's
    # Authorization value with <AUTHORIZATION>
    c.filter_sensitive_data('<AUTHORIZATION>') do |interaction|
      interaction.request.headers['Authorization'].first
    end

    c.before_record do |interaction|
      key = prepare_vcr_secret

      if interaction.response.body && !interaction.response.body.empty?
        cipher = OpenSSL::Cipher.new("AES-256-CBC")
        iv = cipher.random_iv

        cipher.encrypt
        cipher.key = key
        cipher.iv = iv

        encrypted = cipher.update(interaction.response.body)
        encrypted << cipher.final

        interaction.response.body = [iv, encrypted].pack('mm')
      end
    end

    c.before_playback do |interaction|
      key = prepare_vcr_secret

      if interaction.response.body && !interaction.response.body.empty?
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

  require 'multipartable'
  def Multipartable.secure_boundary
    "--THIS-IS-STABLE-FOR-TESTING"
  end
end

def prepare_vcr_secret
  p ENV
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