module Etna
  class Route
    attr_reader :name

    def initialize(method, options, &block)
      @method = method
      if options.is_a?(Hash)
        path, @action = options.first
        @path = path_regexp(path)
        @name = options[:as]
      else
        @path = path_regexp(options)
        @block = block
      end
    end

    def matches? request
      matches_path?(request.path) && @method == request.request_method
    end

    def call(app, request)
      update_params(request)

      if @action
        controller, action = @action.split('#')
        controller_class = Kernel.const_get(
          :"#{controller.camel_case}Controller"
        )
        logger = request.env['rack.logger']
        logger.warn("Calling #{controller}##{action}")
        return controller_class.new(request, action).response
      elsif @block
        return app.instance_eval(&@block)
      end
    end

    def update_params request
      match = @path.match(request.path)
      request.env['rack.request.params'].update(
        Hash[
          match.names.map(&:to_sym).zip(
            match.captures.map do |capture|
              URI.decode(capture)
            end
          )
        ]
      )
    end

    def matches_path? path
      path.match(@path)
    end

    def path_regexp(path)
      Regexp.new(
        path
          .gsub(/:([\w]+)/, '(?<\1>\w+)')
          .gsub(/\*([\w]+)$/, '(?<\1>.*)')
      )
    end
  end
end
