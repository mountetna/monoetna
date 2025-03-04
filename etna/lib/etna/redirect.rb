module Etna
  def self.Redirect(req)
    Redirect.new(req)
  end

  class Redirect
    def initialize(request)
      @request = request
    end

    def to(path, &block)
      return Rack::Response.new(
        [ { errors: [ 'Cannot redirect out of domain' ] }.to_json ], 422,
        { 'Content-Type' => 'application/json' }
      ).finish unless matches_domain?(path)

      response = Rack::Response.new

      response.redirect(path.gsub("http://", "https://"), 302)

      yield response if block_given?

      response.finish
    end

    private

    def matches_domain?(path)
      top_domain(Etna::Application.instance.host) == top_domain(URI(path).host)
    end

    def top_domain(host_name)
      host_name.split('.')[-2..-1].join('.')
    end
  end
end
