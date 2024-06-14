require "connection_pool"

class Vulcan

  class SSHConnectionPool
    def initialize(host, username, password, pool_size: 5, timeout: 300)
      @pool = ConnectionPool.new(size: pool_size, timeout: timeout) do
        Net::SSH.start(host, username, password: password)
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