from airflow.plugins_manager import AirflowPlugin
from etna.plugins.impersonate_user_view import impersonate_user_view


impersonate_package = {
    "name": "Impersonate",
    "category": "Admin",
    "view": impersonate_user_view,
}

class ImpersonateUserPlugin(AirflowPlugin):
    name = 'impersonate_user_plugin'

    appbuilder_views = [impersonate_package] if impersonate_user_view is not None else []