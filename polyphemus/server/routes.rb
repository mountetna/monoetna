Polyphemus = Polyphemus.new()

# Page view end points
Polyphemus.add_route('GET', '/', 'ClientController#basic_view')
Polyphemus.add_route('GET', '/logged-out', 'ClientController#logged_out')
Polyphemus.add_route('GET', '/user-admin', 'ClientController#user_admin')
Polyphemus.add_route('GET', '/network-utils', 'ClientController#network_utils')

# User logging end points
Polyphemus.add_route('POST', '/login', 'UserLogController#log_in')
Polyphemus.add_route('POST', '/logout', 'UserLogController#log_out')
Polyphemus.add_route('POST', '/check', 'UserLogController#check_log')

# User administration end points
Polyphemus.add_route('POST', '/get-users', 'UserAdminController#get_users')
Polyphemus.add_route('POST', '/get-projects', 'UserAdminController#get_projects')
Polyphemus.add_route('POST', '/get-permissions', 'UserAdminController#get_permissions')
Polyphemus.add_route('POST', '/upload-permissions', 'UserAdminController#upload_permissions')
Polyphemus.add_route('POST', '/remove-permissions', 'UserAdminController#remove_permissions')

# Invalidate all tokens in the system.
# This could cause problems with Metis if there are active uploads.
Polyphemus.add_route('POST', '/logout-all', 'UserAdminController#logout_all')