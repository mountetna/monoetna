require "connection_pool"

class Vulcan

  class SSHConnectionPool
    def initialize(host, username, password, use_ssh_config, settings: {}, pool_size: 5, timeout: 600, max_retries: 3, retry_delay: 5)
      @host = host
      @username = username
      @password = password
      @use_ssh_config = use_ssh_config
      @settings = settings.merge(
        keepalive: true,
        keepalive_interval: 30,
        keepalive_maxcount: 3,
        timeout: timeout
      )
      @max_retries = max_retries
      @retry_delay = retry_delay
      
      @pool = ConnectionPool.new(size: pool_size, timeout: timeout) do
        create_connection
      end
    end

    def create_connection
      retries = 0
      begin
        if @password
          Vulcan.instance.logger.info("Establishing SSH connection to #{@host} with username: #{@username} and password ...")
          if @use_ssh_config
            Vulcan.instance.logger.info("Reading ssh config from ~/.ssh/config ...")
            Net::SSH.start(@host, @username, password: @password, config: true, **@settings)
          else
            Net::SSH.start(@host, @username, password: @password, config: false, **@settings)
          end
        else
          Vulcan.instance.logger.info("Establishing SSH connection to #{@host} with user #{@username} and SSH config ...")
          Vulcan.instance.logger.info("Reading SSH config from ~/.ssh/config ...")
          Net::SSH.start(@host, @username, config: true, **@settings)
        end
      rescue Net::SSH::Exception, Errno::ECONNREFUSED, IOError => e
        retries += 1
        if retries <= @max_retries
          Vulcan.instance.logger.warn("SSH connection attempt #{retries}/#{@max_retries} failed: #{e.message}. Retrying in #{@retry_delay} seconds...")
          sleep @retry_delay
          retry
        else
          Vulcan.instance.logger.error("Failed to establish SSH connection after #{@max_retries} attempts: #{e.message}")
          raise e
        end
      end
    end

    def with_conn
      @pool.with do |ssh|
        begin
          begin
            channel = ssh.open_channel do |ch|
              ch.exec("echo 'connection test'") do |ch, success|
                raise "Could not execute command" unless success
              end
            end
            channel.wait
          rescue Net::SSH::Exception, IOError => e
            Vulcan.instance.logger.warn("SSH connection validation failed: #{e.message}")
            raise e
          end

          yield ssh
        rescue Net::SSH::Disconnect, Errno::ECONNREFUSED, IOError => e
          Vulcan.instance.logger.warn("SSH connection lost: #{e.message}, attempting to reconnect...")
          # Force pool to create a new connection
          @pool.shutdown { |conn| conn.close rescue nil }
          ssh = create_connection
          yield ssh
        end
      end
    end

    def close_all
      @pool.shutdown { |conn| conn.close rescue nil }
    end
  end

end