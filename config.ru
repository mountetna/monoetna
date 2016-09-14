# config.ru

require 'rack'
require 'json'
require 'digest'

require './server/conf'
require './server/controller'
require './server/utils'
require './server/metis'
require './server/routes'

use Rack::Static, urls: ['/css', '/js', '/fonts'], root: 'client'

run Metis