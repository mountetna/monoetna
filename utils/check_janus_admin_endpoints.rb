  # This will check the Janus URL endpoints.
require 'net/http'
require 'openssl'
require 'json'

require './conf.rb'
require './utils.rb'
require '../server/secrets'

class JanusCheckAdmin

# Setup a user
  def initialize()

    @good_email = 'jason.cater@ucsf.edu'
    @good_passwd = 'w4yn0w4y'
    @good_token = nil
  end

# Log in the user
  def log_in() 

    data = {

      :email=> @good_email, 
      :pass=> @good_passwd, 
      :app_key=> Secrets::APP_KEY
    }
    url = Conf::JANUS_ADDR+'login'
    response = JSON.parse(Utils::make_request(url, data))
    @good_token = Utils::verify_login(response)
    puts response
    puts ''
  end

# Check the end points
  def get_users()

    pre_request('get-users')
  end

  def get_projects()

    pre_request('get-projects')
  end

  def get_groups()

    pre_request('get-groups')
  end

  def get_permissions()

    pre_request('get-permissions')
  end

  def pre_request(uri)

    url = Conf::JANUS_ADDR+uri

    puts url

    data = { :token=> @good_token, :app_key=> Secrets::APP_KEY }
    response = JSON.parse(Utils::make_request(url, data))
    puts JSON.pretty_generate(response)
    puts ''
  end
end

janusCheckAdmin = JanusCheckAdmin.new()
janusCheckAdmin.log_in()
janusCheckAdmin.get_groups()