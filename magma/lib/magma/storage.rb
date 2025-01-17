class Magma
  class Storage
    def self.setup
      config = Magma.instance.config(:storage)

      return Magma::Storage.new unless config

      case config[:provider]
      when 'aws'
        return Magma::Storage::Aws.new(config)
      when 'metis'
        return Magma::Storage::Metis.new(config)
      end
    end

    def download_url project_name, path
    end

    def upload_url(project_name, path, params={})
    end

    class Metis < Magma::Storage
      def initialize(config)
        @config = config
      end

      def download_url project_name, path
        hmac_url(
          'GET',
          @config[:host],
          "/#{project_name}/download/magma/#{path}",
          @config[:download_expiration]
        )
      end

      def upload_url(project_name, path, params = {})
        hmac_url(
          'POST',
          @config[:host],
          "/#{project_name}/upload/magma/#{path}",
          @config[:upload_expiration],
          params
        )
      end

      private

      def hmac_url(method, host, path, expiration=0, headers={})
        URI::HTTPS.build(
          Etna::Hmac.new(
            Magma.instance,
            method: method,
            host: host,
            path: encode_path(path),
            expiration: (Time.now + expiration).iso8601,
            nonce: Magma.instance.sign.uid,
            id: :magma,
            headers: headers
          ).url_params
        )
      end

      def encode_path(path)
        path.split('/').map { |c| Etna::Route.encode_path_component(c) }.join('/')
      end
    end

    class Aws < Magma::Storage
      def initialize(config)
        require 'fog/aws'
        require 'carrierwave/sequel'
        require 'carrierwave/storage/fog'

        @config = config
        @fog = Fog::Storage.new(@config[:credentials])

        carrier_wave_init
      end

      def download_url project_name, path
        @fog.get_object_url(
          @config[:directory],
          "uploads/#{path}",
          Time.now + @config[:expiration]*60,
          path_style: true
        )
      end

      def upload_url project_name, path
        @fog.put_object_url(
          @config[:directory],
          "uploads/#{path}",
          Time.now + @config[:expiration]*60,
          path_style: true
        )
      end

      private

      def carrier_wave_init
        CarrierWave.tmp_path = Magma.instance.config(:tmp_path)
        CarrierWave.configure do |config|
          config.storage :fog
          config.fog_provider = 'fog/aws'
          config.fog_credentials = @config[:credentials]
          config.fog_directory = @config[:directory]
          config.fog_public = false
          config.fog_attributes = {'Cache-Control'=> "max-age=#{365 * 86400}"}
        end
      end
    end
  end
end
