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
      @match_ext = options[:match_ext]
    end

    def to_hash
      {
        method: @method,
        route: @route,
        name: @name.to_s,
        params: parts
      }.compact
    end

    def matches?(request)
      @method == request.request_method && request.path.match(route_regexp)
    end

    NAMED_PARAM=/:([\w]+)/
    GLOB_PARAM=/\*([\w]+)$/

    PARAM_TYPES=[ NAMED_PARAM, GLOB_PARAM ]

    UNSAFE=/[^\-_.!~*'()a-zA-Z\d;\/?:@&=+$,]/

    def self.path(route, params=nil)
      if params
        PARAM_TYPES.reduce(route) do |path,pat|
          path.gsub(pat) do
           URI.encode( params[$1.to_sym], UNSAFE)
          end
        end
      else
        route
      end
    end

    def path(params=nil)
      self.class.path(@route, params)
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
        return [ 403, { 'Content-Type' => 'application/json' }, [ { error: 'You are forbidden from performing this action.' }.to_json ] ]
      end

      if @action
        controller, action = @action.split('#')
        controller_class = Kernel.const_get(
          :"#{controller.camel_case}Controller"
        )
        logger = request.env['etna.logger']
        user = request.env['etna.user']

        params = request.env['rack.request.params'].map do |key,value|
          [ key, redact(key, value) ]
        end.to_h

        logger.warn("User #{user ? user.email : :unknown} calling #{controller}##{action} with params #{params}")
        return controller_class.new(request, action).response
      elsif @block
        application = Etna::Application.find(app.class).class
        controller_class = application.const_defined?(:Controller) ? application.const_get(:Controller) : Etna::Controller

        controller_class.new(request).response(&@block)
      end
    end

    # the route does not require authorization
    def noauth?
      @auth && @auth[:noauth]
    end

    def ignore_janus?
      @auth && @auth[:ignore_janus]
    end

    private

    def application
      @application ||= Etna::Application.instance
    end

    def compact(value)
      value = value.to_s
      value = value[0..500] + "..." + value[-100..-1] if value.length > 600
      value
    end

    def redact_keys
      @redact_keys ||= application.config(:log_redact_keys).split(",").map do |key|
        key.to_sym
      end
    end

    def redact(key, value)
      # From configuration, redact any values for the supplied key values, so they
      #   don't appear in the logs.
      return compact(value) unless application.config(:log_redact_keys)

      if value.is_a?(Hash)
        redacted_value = value.map do |value_key, value_value|
          [ value_key, redact(value_key, value_value) ]
        end.to_h
        return redacted_value
      elsif redact_keys.include?(key)
        return "*"
      end

      return compact(value)
    end

    def authorized?(request)
      # If there is no @auth requirement, they are ok - this doesn't preclude
      # them being rejected in the controller response
      !@auth || @auth[:noauth] || (user_authorized?(request) && hmac_authorized?(request))
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
      # either there is no hmac requirement, or we have a valid hmac
      !@auth[:hmac] || request.env['etna.hmac']&.valid?
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

    def separator_free_match
      if @match_ext
        '(?<\1>[^\/\?]+)'
      else
        '(?<\1>[^\.\/\?]+)'
      end
    end

    def route_regexp
      @route_regexp ||=
        Regexp.new(
          '\A' +
          @route.
            # any :params match separator-free strings
            gsub(NAMED_PARAM, separator_free_match).
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
