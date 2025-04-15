require "connection_pool"

class Vulcan

  class SSHConnectionPool
    def initialize(host, username, password, use_ssh_config, settings: {}, pool_size: 5, timeout: 300)
      @pool = ConnectionPool.new(size: pool_size, timeout: timeout) do
        if password
          Vulcan.instance.logger.info("Establishing SSH connection to #{host} with username: #{username} and password ...")
          # Direct connection using username and password
          if use_ssh_config
            Vulcan.instance.logger.info("Reading ssh config from ~/.ssh/config ...")
            Net::SSH.start(host, username, password: password, config: true, **settings)
          else
            Net::SSH.start(host, username, password: password, config: false, **settings)
          end
        else
          Vulcan.instance.logger.info("Establishing SSH connection to #{host} with user #{username} and SSH config ...")
          Vulcan.instance.logger.info("Reading SSH config from ~/.ssh/config ...")
          # Use SSH config for connection
          Net::SSH.start(host, username, config: true, **settings)
        end
      end
    end

    def with_conn
      @pool.with do |ssh|
        yield ssh
      end
    end

    def close_all
      @pool.shutdown { |conn| conn.close }
    end
  end

end