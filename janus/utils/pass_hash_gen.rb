#!/home/developer/.rbenv/shims/ruby

require 'digest'
require 'yaml'
require '../lib/server/service/sign_service'

def print_help
  puts ''
  puts ' desc: A simple util to generate a passhash for the Janus Auth server.'
  puts 'usage: ./pass_hash_gen.rb [environment] [plain text password]'
  puts ''
end

# If the argument is 'nil', or if there are to many arguments, or if the
# argument is an empty string, then bail.
if ARGV[0].nil? || ARGV.length != 2 || ARGV[0].to_s.length == 0
  print_help
  exit 1
end

config = YAML::load(File.open('../config.yml'))
if !config.key?(ARGV[0].to_sym)
  puts 'ERROR: Specify an environment to generate a key against. The following'\
' enviroments are defined in config.yml:'
  puts config.keys
  exit 1
end

puts SignService::hash_password(
  [
    ARGV[1].to_s,
    config[ARGV[0].to_sym][:pass_salt]
  ],
  config[ARGV[0].to_sym][:pass_algo]
)