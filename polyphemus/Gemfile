source 'https://rubygems.org'

gem 'rack'
gem 'sequel'
gem 'pg'
gem 'puma', '5.0.2'
gem 'etna', git: 'https://github.com/mountetna/monoetna.git', branch: 'refs/artifacts/gem-etna/60667019a69dda8b6246a241d9962b717081b9bd'
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
