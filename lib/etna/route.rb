module Etna
  class Route
    def initialize(path, method, action, &block)
      @path = path
      @method = method
      @action = action
      @block = block
    end

    def matches? request
      @path == request.path && @method == request.request_method
    end

    def call(app, request)
      if @action
        controller, action = @action.split('#')
        controller_class = Kernel.const_get(
          :"#{controller.camel_case}Controller"
        )
        logger = request.env['rack.logger']
        logger.warn("Calling #{controller}##{action}")
        return controller_class.new(request, action).response
      elsif @block
        return app.instance_eval(@block)
      end
    end
  end
end
