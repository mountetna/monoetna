require('rack')
require('json')
require('pg')
require('sequel')
require('digest')
require('fileutils')
require('net/http')
require('logger')
require('openssl')
require('securerandom') # Required by Phusion Passenger

# The details
require('./server/conf')
require('./server/secrets')

# The database
require('./server/service/postgres_service')
PostgresService::connect
require('./server/models/file_model')

# The application
require('./server/models/user_model')
require('./server/service/sign_service')
require('./server/errors/basic_error')
require('./server/metis')
require('./server/routes')
require('./server/controllers/basic_controller')
require('./server/controllers/main_controller')
require('./server/controllers/upload_controller')
require('./server/controllers/client_controller')
require('./server/controllers/user_log_controller')

use(Rack::Static, {urls: ['/css', '/js', '/fonts', '/img'], root: 'client'})
run(Metis)
