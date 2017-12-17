module Etna
  class Route
    attr_reader :name

    def initialize(method, route, options, &block)
      @method = method
      @action = options[:action]
      @name = options[:as]
      @route = route.gsub(/\A(?=[^\/])/, '/')
      @block = block
    end

    def matches?(request)
      @method == request.request_method && request.path.match(route_regexp)
    end

    def path(params)
      @route
        .gsub(/:([\w]+)/) { params[$1.to_sym] }
        .gsub(/\*([\w]+)$/) { params[$1.to_sym] }
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

    private

    def update_params(request)
      match = route_regexp.match(request.path)
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

    def route_regexp
      @route_regexp ||=
        Regexp.new(
          '\A' +
          @route.
            # any :params match separator-free strings
            gsub(/:([\w]+)/, '(?<\1>[^\.\/\?]+)').
            # any *params match arbitrary strings
            gsub(/\*([\w]+)$/, '(?<\1>.+)').
            # ignore any trailing slashes in the route
            gsub(/\/\z/, '') +
          # trailing slashes in the path can be ignored
          '/?' +
          '\z'
        )
    end
  end
end
