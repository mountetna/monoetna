# redis_service.rb

class RedisService

  def initialize()

    @redis = Redis.new(:host => 'localhost', :port => 6379, :db => 0)
    
    if @redis.get('index') == nil

      # Replace this with a bit of code that finds the maximum index value in
      # the keys.
      @redis.set('index', 0)
    end 
  end

  def get_new_index()

    @redis.incr('index')
    index = @redis.get('index')
    index = index.to_i.to_s(16)
    
    if(index.length < 8)

      pad_length = 8 - index.length
      index = ('0' * pad_length) + index
    end

    return index
  end

  def retrieve_file_key(pattern)

    @redis.keys(pattern)
  end

  def retrieve_file_status(key)

    file_status = @redis.get(key)
    return (file_status == nil) ? nil : JSON.parse(file_status)
  end

  def set_file_status(key, value)

    @redis.set(key, value)
  end

  def remove_file_status(key)

    @redis.del(key)
  end

  def status_present?(key)

    return (@redis.keys(key).length == 0) ? false : true
  end
end