![Run JavaScript tests](https://github.com/mountetna/etna/workflows/Run%20JavaScript%20tests/badge.svg)  ![Run Ruby tests](https://github.com/mountetna/etna/workflows/Run%20Ruby%20tests/badge.svg)



# Mount Etna common library

This repository is meant to provide a Ruby gem and some JavaScript
modules to give a common server infrastructure for Mount Etna projects.

Documentation for the Etna framework can be found at https://mountetna.github.io/

## Ruby Gem

The gem may be built by hand or installed from github with bundler:

    gem 'etna', github: 'mountetna/etna'

A basic Mount Etna application uses a Rack server and a Sequel database with a
postgresql adapter.



