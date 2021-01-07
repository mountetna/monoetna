source 'https://rubygems.org'

gem 'rack'
gem 'sequel'
gem 'pg'
gem 'nokogiri'
gem 'puma', '5.0.2'
gem 'etna', git: 'https://github.com/mountetna/monoetna.git', branch: 'refs/artifacts/gem-etna/f81c8215a903a130a10826b9904cf9d8f7b7d4d0'
gem 'actionpack' # For streaming the job controller results back...

group :development, :test do
  gem 'rack-test', require: "rack/test"
  gem 'rspec'
  gem 'simplecov'
  gem 'pry'
  gem 'pry-byebug'
  gem 'webmock'
  gem 'debase'
  gem 'database_cleaner'
  gem 'vcr'
end
