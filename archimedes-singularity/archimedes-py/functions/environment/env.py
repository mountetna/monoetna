import os

token = os.getenv('TOKEN', None)
# In future, convert from the user's standard token to a longer-lived but read only, and this project only, version.

metis_host = os.getenv('METIS_HOST', None)
magma_host = os.getenv('MAGMA_HOST', None)
project_name = os.getenv('PROJECT_NAME', None)
