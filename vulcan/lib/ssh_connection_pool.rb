require "connection_pool"
require "net/ssh"

class Vulcan
  class SSHConnectionPool
    def initialize(host, username, password, use_ssh_config,
                   settings: {}, pool_size: 5, timeout: 600,
                   max_retries: 3, retry_delay: 5)

      @host           = host
      @username       = username
      @password       = password
      @use_ssh_config = use_ssh_config

      # 1) Tune keep-alive for faster dead-socket detection (~30s total)
      @settings = settings.merge(
        keepalive:          true,
        keepalive_interval: 10,
        keepalive_maxcount: 3,
        timeout:            timeout
      )

      @max_retries  = max_retries
      @retry_delay  = retry_delay
      @pool_size    = pool_size
      @pool_timeout = timeout

      # 2) Build initial pool
      @pool = build_pool
    end

    # Borrow a session and yield; on session-level errors, rebuild and retry a few times.
    def with_conn
      attempts = 0
      begin
        @pool.with do |ssh|
          yield ssh
        end
      rescue Net::SSH::AuthenticationFailed, Net::SSH::HostKeyMismatch => e
        Vulcan.instance.logger.error("SSH auth/hostkey error: #{e.class}: #{e.message}")
        raise
      rescue Net::SSH::Disconnect,Errno::ECONNREFUSED, Errno::ECONNRESET, Errno::EBADF, IOError, SocketError, Timeout::Error => e
        attempts += 1
        if attempts <= @max_retries
          Vulcan.instance.logger.warn(
            "SSH session error (#{e.class}: #{e.message}), attempt #{attempts}/#{@max_retries}. "\
            "Recreating pool and retrying in #{@retry_delay}s…"
          )
          recreate_pool_once
          sleep @retry_delay
          retry
        else
          Vulcan.instance.logger.error("Giving up after #{attempts} attempts: #{e.class}: #{e.message}")
          raise
        end
      end
    end

    def recreate_pool_once
      @rebuild_mutex ||= Mutex.new
      @rebuild_mutex.synchronize do
        Vulcan.instance.logger.info("Recreating SSH connection pool for #{@host}…")
        @pool.shutdown { |conn| conn.close rescue nil }
        @pool = build_pool
      end
    end

    def close_all
      @pool.shutdown { |conn| conn.close rescue nil }
    end

    private

    def build_pool
      ConnectionPool.new(size: @pool_size, timeout: @pool_timeout) do
        create_connection
      end
    end


    def create_connection
      tries = 0
      begin
        Vulcan.instance.logger.info("Establishing SSH to #{@host} as #{@username}…")
        opts = @settings.dup
        opts[:config] = @use_ssh_config
        opts[:password] = @password if @password
        Net::SSH.start(@host, @username, **opts)
      rescue Net::SSH::Exception, Errno::ECONNREFUSED, Errno::ECONNRESET, IOError, Errno::ETIMEDOUT, SocketError => e
        tries += 1
        if tries <= @max_retries
          Vulcan.instance.logger.warn(
            "SSH connect attempt #{tries}/#{@max_retries} failed: #{e.message}. "\
            "Retrying in #{@retry_delay}s…"
          )
          sleep @retry_delay
          retry
        else
          Vulcan.instance.logger.error(
            "Failed to establish SSH after #{@max_retries} attempts: #{e.message}"
          )
          raise
        end
      end
    end
  end
end
