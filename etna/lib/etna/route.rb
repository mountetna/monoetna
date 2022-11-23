require 'digest'
require 'date'
require "addressable/uri"
require_relative "./censor"
require_relative './instrumentation'

module Etna
  class Route
    include Etna::Instrumentation

    attr_reader :name

    def initialize(method, route, options, &block)
      @method = method
      @action = options[:action]
      @auth = options[:auth]
      @name = route_name(options)
      @route = route.gsub(/\A(?=[^\/])/, '/')
      @block = block
      @match_ext = options[:match_ext]
      @log_redact_keys = options[:log_redact_keys]
      @dont_log = options[:dont_log]
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

    def self.path(route, params=nil, &block)
      if params
        PARAM_TYPES.reduce(route) do |path,pat|
          path.gsub(pat) do
            params[$1.to_sym].split('/').map { |c| block_given? ? yield(c) : Addressable::URI.normalized_encode(c) }.join('/')
          end
        end
      else
        route
      end
    end

    def path(params=nil, &block)
      self.class.path(@route, params, &block)
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

      tags = {}
      if @action
        tags[:controller], tags[:action] = @action.split('#')
      elsif @name
        tags[:action] = @name
      else
        tags[:action] = @route
      end

      params = request.env['rack.request.params']
      user = request.env['etna.user']
      tags[:user_hash] = user ? hash_user_email(user.email) : 'unknown'
      if params && (params.include?(:project_name) || params.include?('project_name'))
        tags[:project_name] = params[:project_name] || params['project_name']
      end

      with_yabeda_tags(tags) do
        increment_it { Yabeda.etna.visits }
        process_call(app, request)
      end
    end

    def hash_user_email(email)
      secret = Etna::Application.instance.config(:user_hash_secret) || 'notsosecret'
      digest = email + secret + Date.today.to_s

      if @name
        digest += @name.to_s
      else
        digest += @route.to_s
      end

      Digest::MD5.hexdigest(digest)
    end

    def process_call(app, request)
      unless authorized?(request)
        if cc_available?(request)
          if request.content_type == 'application/json'
            return [ 403, { 'Content-Type' => 'application/json' }, [ { error: 'You are forbidden from performing this action, but you can visit the project home page and request access.' }.to_json ] ]
          else
            return cc_redirect(request)
          end
        end

        return [ 403, { 'Content-Type' => 'application/json' }, [ { error: 'You are forbidden from performing this action.' }.to_json ] ]
      end

      request.env['etna.redact_keys'] = @log_redact_keys
      request.env['etna.dont_log'] = @dont_log

      if @action
        controller, action = @action.split('#')
        controller_class = Kernel.const_get(
          :"#{controller.camel_case}Controller"
        )
        logger = request.env['etna.logger']
        user = request.env['etna.user']

        return controller_class.new(request, action).response
      elsif @block
        application = Etna::Application.find(app.class).class
        controller_class = application.const_defined?(:Controller) ? application.const_get(:Controller) : Etna::Controller

        controller_class.new(request).response(&@block)
      end
    end

    time_it(:process_call) { Yabeda.etna.response_time }

    # the route does not require authorization
    def noauth?
      @auth && @auth[:noauth]
    end

    def ignore_janus?
      @auth && @auth[:ignore_janus]
    end

    def has_user_constraint?(constraint)
      @auth && @auth[:user] && @auth[:user][constraint]
    end

    private

    def janus
      @janus ||= Etna::JanusUtils.new
    end

    # If the application asks for a code of conduct redirect
    def cc_redirect(request, msg = 'You are unauthorized')
      return [ 401, { 'Content-Type' => 'text/html' }, [msg] ] unless application.config(:auth_redirect)

      params = request.env['rack.request.params']

      uri = URI(
        application.config(:auth_redirect).chomp('/') + "/#{params[:project_name]}/cc"
      )
      uri.query = URI.encode_www_form(refer: request.url)
      Etna::Redirect(request).to(uri.to_s)
    end

    def application
      @application ||= Etna::Application.instance
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
        user.respond_to?(constraint) && (
          param_name.is_a?(Symbol) ?
            user.send(constraint, params[param_name]) :
            user.send(constraint, param_name))
      end
    end

    def cc_available?(request)
      user = request.env['etna.user']

      return false unless user

      params = request.env['rack.request.params']

      return false unless params[:project_name]

      # Only check for a CC if the user does not currently have permissions
      #   for the project
      return false if user.permissions.keys.include?(params[:project_name])

      !janus.community_projects(user.token).select do |project|
        project.project_name == params[:project_name]
      end.first.nil?
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
              capture.split('/').map {|c| Addressable::URI.unencode_component(c) }.join('/')
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
