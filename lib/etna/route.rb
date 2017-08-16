module Etna
  class Route
    def initialize method, action, &block
      @method = method
      @action = action
      @block = block
    end

    def call(app, request)
      if @action
        controller, action = @action.split('#')
        controller_class = Kernel.const_get(
          :"#{controller.camel_case}Controller"
        )
        app.logger.warn("Calling #{controller}##{action}")
        return controller_class.new(request, action).response
      elsif @block
        return app.instance_eval(@block)
      end
    end
  end
end
