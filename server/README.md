## Some notes about setting up.

You are going to need a `secrets.rb` file which will contain your app secrets 
and it should look like so...

./metis/server/secrets.rb

```
module Secrets

  SECRET_KEY = [A KEY USED FOR HMAC AUTHENTICATIONS]
  APP_KEY = [THE APP AUTH KEY IN THE JANUS AUTH SERVER]

  PSQL_USER = [POSTGRES USER]
  PSQL_PASS = [POSTGRES PASSWORD]

  # This is the root directory for the data.
  ROOT_DIR = [THE DIR WHERE THE FILES ARE KEPT]

  # This is the address for the Janus auth server.
  JANUS_ADDR = [ADDR OF THE AUTH SERVER]
end
```