# Policyfile.rb - Describe how you want Chef Infra Client to build your system.
#
# For more information on the Policyfile feature, visit
# https://docs.chef.io/policyfile.html

# Specify a custom source for a single cookbook:
cookbook 'mountetna', path: '.'
cookbook 'docker', path: '../docker'
