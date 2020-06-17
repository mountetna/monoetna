## Some notes about setting up.

You are going to need a `secrets.rb` file which will contain your app secrets 
and it should look like so...

./polyphemus/server/secrets.rb


```
module Secrets

  APP_KEY = '[Polyphemus App Key]'
  JANUS_ADDR = 'http://[Janus URL]'
end
```