module Etna
  class Command
    class << self
      def usage(desc)
        define_method :usage do
          "  #{"%-30s" % name}#{desc}"
        end
      end
    end

    def name
      self.class.name.snake_case.split(/::/).last.to_sym
    end

    # To be overridden during inheritance.
    def execute
      raise "Command is not implemented"
    end

    # To be overridden during inheritance, to e.g.
    # connect to a database.
    # Should be called with super by inheriting method
    def setup(config)
      Etna::Application.find(self.class).configure(config)
    end
  end
end
