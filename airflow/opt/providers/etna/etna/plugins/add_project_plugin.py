from airflow.plugins_manager import AirflowPlugin
from etna.plugins.add_project_view import add_project_view


add_project_package = {
    "name": "Add Project",
    "category": "Admin",
    "view": add_project_view,
}

class AddProjectPlugin(AirflowPlugin):
    name = 'add_project_plugin'

    appbuilder_views = [add_project_package]
