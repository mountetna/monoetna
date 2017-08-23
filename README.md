# Mount Etna common library

This repository is meant to provide a Ruby gem (and later some javascript
modules) to give a common server infrastructure for Mount Etna projects.

The gem may be built by hand or installed from github with bundler:

    gem 'etna', github: 'mountetna/etna'

A basic Mount Etna application uses a Rack server and a Sequel database with a
postgresql adapter.

## Etna::Application

The base class of your application (e.g. 'Ariadne') should include Etna::Application

    class Ariadne
      include Etna::Application
    end

This will make the class a singleton (using Ruby's Singleton interface). The main purpose
of Etna::Application is to hold configuration and make a connection to a database, separate
from the Rack server (so you can, e.g., connect to the database via a command-line interface).
Configuration is a Hash object usually stored in a config.yml at the root of the application:

    Ariadne.instance.configure(YAML.load(File.read("config.yml")))

Configurations are separated out by environment, e.g. your config.yml might look like this:

    :test:
      :log_path: ./test_log

    :development:
      :log_path: ./log

The current environment is set using an environment variable based on your application name. (e.g. ENV["ARIADNE_ENV"]). You may access configuration for the current environment from your singleton:

    Ariadne.instance.config(:log_path)

##Etna::Server

This class implements a basic Rack server and provides simple routing, logging, and error-handling.

The server should use the application class as a container:

    class Ariadne
      class Server < Etna::Server
      end
    end

Routes may be either defined using 'controller#action' syntax or using a block:

    class Ariadne
      class Server < Etna::Server
        route '/weave', 'loom#weave'

	route '/trash-talk' do
	  TrashTalkController.new(@request).response
	end
      end
    end

###Etna::Controller

Controllers are in the global namespace and should be named like SomethingController. Actions are
methods defined on the controller.

    class LoomController < Etna::Controller
      def weave
        success('application/json', tapestry: 'image-of-aphrodite.jpg')
      end
    end

The controller provides @request and @response rack objects, @params (built by Etna::ParseBody), and #log, #success and #failure methods.

The action is usually invoked by the #response method. Both the action and ultimately the #response method must return a valid Rack response (e.g. [ 200, {}, 'OK' ])

###Etna::Error

There are two error classes available: Etna::BadRequest (for client errors) and Etna::ServerError (for server errors).
Raising will cause the controller to return a status of 422 or 500 with the given error.
