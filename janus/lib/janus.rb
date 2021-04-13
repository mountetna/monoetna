require 'sequel'
require 'active_support/all'

class Janus
  include Etna::Application
  attr_reader :db

  def setup_db(load_models=true)
    @db = Sequel.connect(config(:db))
    @db.extension :connection_validator
    @db.extension :pg_json
    @db.pool.connection_validation_timeout = -1

    require_relative 'models' if load_models
  end

  def set_token_cookie(response,token)
    # Tear apart token to get expire time
    expire_time = Time.at(
      JSON.parse(
        Base64.decode64(
          token.split('.')[1]
        )
      )["exp"]
    )

    # Set cookie
    response.set_cookie(
      config(:token_name),
      value: token,
      path: '/',
      domain: config(:token_domain),
      expires: expire_time,
      secure: true,
      same_site: :strict
    )
  end

end
