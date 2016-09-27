# redis_service.rb

class RedisService

  def initialize()

    @redis = Redis.new(:host => 'localhost', :port => 6379, :db => 0)
  end

  def retrive_file_status(key)
    
    @redis.get(key)
  end

  def set_file_status(key, value)

    @redis.set(key, value)
  end

  def status_present?(key)

    if @redis.keys(key).length == 0

      return false
    else

      return true
    end
  end
end