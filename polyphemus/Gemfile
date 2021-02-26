source 'https://rubygems.org'

gem 'rack'
gem 'sequel'
gem 'pg'
gem 'puma', '5.0.2'
gem 'etna', git: 'https://github.com/mountetna/monoetna.git', branch: 'refs/artifacts/gem-etna/e95581f0ec78dff9f9ab6997bab6212358cfe9cf'
gem 'actionpack' # For streaming the job controller results back...
gem 'aspera-cli'
gem 'concurrent-ruby'

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
