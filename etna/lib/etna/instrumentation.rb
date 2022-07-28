module Etna
  module Instrumentation
    def self.included(cls)
      cls.instance_eval do
        def self.time_it(method_name, &metric_block)
          orig_method_name = :"#{method_name}_without_time_it"
          self.alias_method orig_method_name, method_name

          self.define_method method_name do |*args, **kwds|
            time_it(orig_method_name, *args, **kwds, &metric_block)
          end
        end
      end
    end

    def increment_it(&block)
      if has_yabeda?
        metric = yield
        metric.increment({})
      end
    end

    def time_it(method_name, *args, **kwds, &metric_block)
      if has_yabeda?
        start = Process.clock_gettime(Process::CLOCK_MONOTONIC)
        begin
          return send(method_name, *args, **kwds)
        ensure
          dur = Process.clock_gettime(Process::CLOCK_MONOTONIC) - start
          if block_given?
            metric = yield
            tags = {}
          else
            tags = {class_name: self.class.name, method_name: method_name.to_s}
            metric = Yabeda.etna.perf
          end
          metric.measure(tags, dur)
        end
      else
        return send(method_name, *args, **kwds)
      end
    end

    def with_yabeda_tags(tags, &block)
      if has_yabeda?
        Yabeda.with_tags(tags, &block)
      else
        yield
      end
    end

    def has_yabeda?
      defined?(Yabeda) && Yabeda.configured?
    end
  end
end