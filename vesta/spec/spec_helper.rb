require "bundler"
require "rack"
require "rack/test"
Bundler.require(:default, :test)

ENV["VESTA_ENV"] = "test"

require "webmock/rspec"
require "database_cleaner"
require "factory_bot"
require "simplecov"
SimpleCov.start

require "yaml"
require "etna/spec/vcr"

require "fileutils"
require "timecop"

require_relative "../lib/server"
require_relative "../lib/vesta"
