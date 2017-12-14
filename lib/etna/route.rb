module Etna
  class Route
    attr_reader :name, :path

    def initialize(method, path, options, &block)
      @method = method
      @action = options[:action]
      @name = options[:as]
      @path = path.gsub(/\A(?=[^\/])/, '/')
      @block = block
    end

    def matches?(request)
      matches_path?(request.path) && @method == request.request_method
    end

    def call(app, request)
      update_params(request)

      if @action
        controller, action = @action.split('#')
        controller_class = Kernel.const_get(
          :"#{controller.camel_case}Controller"
        )
        logger = request.env['etna.logger']
        logger.warn("Calling #{controller}##{action}")
        return controller_class.new(request, action).response
      elsif @block
        Etna::Controller.new(request).instance_eval(&@block)
      end
    end

    def update_params(request)
      match = path_regexp.match(request.path)
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

    def matches_path?(path)
      path.match(path_regexp)
    end

    def path_regexp
      Regexp.new(
        '\A' +
        @path
          .gsub(/:([\w]+)/, '(?<\1>\w+)')
          .gsub(/\*([\w]+)$/, '(?<\1>.*)') +
        '\z'
      )
    end
  end
end
