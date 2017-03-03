# The packages
require 'rack'
require 'json'
require 'digest'
require 'logger'

# The details
require './server/conf'
require './server/secrets'

# The application
require './server/polyphemus'
require './server/routes'
require './server/controllers/client_controller'
require './server/controllers/network_utils_controller'
require './server/controllers/user_admin_controller'

use Rack::Static, urls: ['/css', '/js', '/fonts', '/img'], root: 'client'
run(Polyphemus)