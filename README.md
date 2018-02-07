# Mount Etna common library

This repository is meant to provide a Ruby gem (and later some javascript
modules) to give a common server infrastructure for Mount Etna projects.

The gem may be built by hand or installed from github with bundler:

    gem 'etna', github: 'mountetna/etna'

A basic Mount Etna application uses a Rack server and a Sequel database with a
postgresql adapter.

## Etna::Application

The base class of your application (e.g. 'Arachne') should include Etna::Application

    class Arachne
      include Etna::Application
    end

This will make the class a singleton (using Ruby's Singleton interface). The main purpose
of Etna::Application is to hold configuration and make a connection to a database, separate
from the Rack server (so you can, e.g., connect to the database via a command-line interface).
Configuration is a Hash object usually stored in a config.yml at the root of the application:

    Arachne.instance.configure(YAML.load(File.read("config.yml")))

Configurations are separated out by environment, e.g. your config.yml might look like this:

    :test:
      :log_path: ./test_log

    :development:
      :log_path: ./log

The current environment is set using an environment variable based on your application name. (e.g. ENV["ARIADNE_ENV"]). You may access configuration for the current environment from your singleton:

    Arachne.instance.config(:log_path)

## Etna::Server

This class implements a basic Rack server and provides simple routing, logging, and error-handling.

The server should use the application class as a container:

    class Arachne
      class Server < Etna::Server
      end
    end

Routes may be either defined using 'controller#action' syntax or using a block:

    class Arachne
      class Server < Etna::Server
        route '/weave', action: 'loom#weave'

        route '/trash-talk' do
          TrashTalk.new(@request).response
        end
      end
    end

### Etna::Controller

Controllers are in the global namespace and should be named like SomethingController. Actions are
methods defined on the controller.

    class LoomController < Etna::Controller
      def weave
        success('application/json', tapestry: 'image-of-aphrodite.jpg')
      end
    end

The controller provides @request and @response rack objects, @params (built by Etna::ParseBody), @server, @logger and @user objects.

The action is usually invoked by the #response method. Both the action and ultimately the #response method must return a valid Rack response (e.g. [ 200, {}, 'OK' ])

The controller also gives several helper methods:
* `view(name)` - returns an HTML view from a file found in the VIEW_PATH
* `erb_view(name)` - returns an HTML view from an ERB template using the controller context.
* `success(content_type, msg)` - a successful rack response
* `failure(status,msg)` - an unsuccessful rack response
* `require_params(*params)` - raise Etna::BadRequest unless @params has params
* `route_path(name, params)` - The URL corresponding to a named route with the given parameters

### Etna::Error

There are two error classes available: Etna::BadRequest (for client errors) and Etna::ServerError (for server errors).
Raising will cause the controller to return a status of 422 or 500 with the given error.

### Etna::ParseBody and Etna::SymbolizeParams

These two Rack layers should be included - the first parses application/json
and multipart messages and makes them available in the rack Request object
under rack.request.params - the Etna::Controller makes this available using the
@params hash. SymbolizeParams turns the keys of the params hash into symbols.

### Etna::Auth

If this middleware is included, it will confirm that the auth token (from an `Authorization: Etna <token>` header or a cookie named in `config(:token_name)`) is a valid JWT verified by the key in `config(:rsa_public)`. If so, it will set 'etna.user' in the Rack environment, which will be picked up by the Etna::Controller as @user. Note that the `:rsa_public` you set in your `config.yml` must correspond to the `:rsa_private` set in your authorization service.

Alternatively, if there is no token, the request is signed by HMAC. HMAC-generated URLs are usually made by the application; there is an Etna::Hmac utility which helps you generate the appropriate requirements and a URL. The result of the HMAC is access to the path with 'etna.hmac' set to `true`.

While you are free to test @user in your controller (e.g., `return failure(403, 'You cannot perform this operation') unless @user.can_edit?(@params[:project_name])`), Etna::Route allows you to do these basic authorization checks via the route. E.g.:

    class Arachne
      class Server < Etna::Server
        get '/test', auth: { user: { can_edit?: :project_name } }
      end
    end

HMAC routes should use `auth: { hmac: true }`.

### Etna::TestAuth

This middleware can be included in lieu of the above to set user credentials without a properly-signed token (useful in testing environments).

### Etna::User

This class validates permissions from the authentication service's JWT. It is not a substitute for a user model, but may be used alongside one.

## Etna::Command

A basic command-line interface can be built with Etna::Command. The Etna::Application singleton will list commands:

    Arachne.instance.commands

and can be setup to run a command from a shell script:

    Arachne.instance.run_command(config, *ARGV)

Commands are subclassed from Etna::Command

    class Arachne
      class Console < Etna::Command
        usage "ariadne console"

        def execute
          require 'irb'
          ARGV.clear
          IRB.start
        end
      end
    end

The command is responsible for creating a database
connection, etc., if appropriate by overriding its #setup
method.
