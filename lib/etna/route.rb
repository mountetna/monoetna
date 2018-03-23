module Etna
  class Route
    attr_reader :name

    def initialize(method, route, options, &block)
      @method = method
      @action = options[:action]
      @auth = options[:auth]
      @name = route_name(options)
      @route = route.gsub(/\A(?=[^\/])/, '/')
      @block = block
    end

    def matches?(request)
      @method == request.request_method && request.path.match(route_regexp)
    end

    NAMED_PARAM=/:([\w]+)/
    GLOB_PARAM=/\*([\w]+)$/

    PARAM_TYPES=[ NAMED_PARAM, GLOB_PARAM ]

    UNSAFE=/[^\-_.!~*'()a-zA-Z\d;\/?:@&=+$,]/

    def path(params=nil)
      if params
        PARAM_TYPES.reduce(@route) do |route,pat|
          route.gsub(pat) do
           URI.encode( params[$1.to_sym], UNSAFE)
          end
        end
      else
        @route
      end
    end

    def parts
      part_list = PARAM_TYPES.map do |pat|
        "(?:#{pat.source})"
      end
      @route.scan(
        /(?:#{part_list.join('|')})/
      ).flatten.compact
    end

    def call(app, request)
      update_params(request)

      unless authorized?(request)
        return [ 403, {}, ['You are forbidden from performing this action.'] ]
      end

      if @action
        controller, action = @action.split('#')
        controller_class = Kernel.const_get(
          :"#{controller.camel_case}Controller"
        )
        logger = request.env['etna.logger']
        logger.warn("Calling #{controller}##{action}")
        return controller_class.new(request, action).response
      elsif @block
        application = Etna::Application.find(app.class).class
        controller_class = application.const_defined?(:Controller) ? application.const_get(:Controller) : Etna::Controller

        controller_class.new(request).response(&@block)
      end
    end

    private

    def authorized?(request)
      # If there is no @auth requirement, they are ok - this doesn't preclude
      # them being rejected in the controller response
      !@auth || (user_authorized?(request) && hmac_authorized?(request))
    end

    def user_authorized?(request)
      # this is true if there are no user requirements
      return true unless @auth[:user]

      user = request.env['etna.user']

      # if there is a user requirement, we must have a user
      return false unless user

      params = request.env['rack.request.params']

      @auth[:user].all? do |constraint, param_name|
        user.respond_to?(constraint) && user.send(constraint, params[param_name])
      end
    end

    def hmac_authorized?(request)
      # either there is no hmac requirement, or we have an hmac
      !@auth[:hmac] || request.env['etna.hmac']
    end

    def route_name(options)
      # use the given one if you can
      return options[:as] if options[:as]

      # otherwise formulate it from the action if possible
      return options[:action].sub(/#/,'_').to_sym if options[:action]

      # unnamed route
      return nil
    end

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
            gsub(NAMED_PARAM, '(?<\1>[^\.\/\?]+)').
            # any *params match arbitrary strings
            gsub(GLOB_PARAM, '(?<\1>.+)').
            # ignore any trailing slashes in the route
            gsub(/\/\z/, '') +
          # trailing slashes in the path can be ignored
          '/?' +
          '\z'
        )
    end
  end
end
