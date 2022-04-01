from airflow.plugins_manager import AirflowPlugin

from etna.hooks.git import GitHook

import logging


# Creating flask appbuilder Menu Items
appbuilder_mitem_toplevel = {
    "name": "Etna Airflow Docs",
    "href": "https://www.apache.org/",
}


class EtnaPlugin(AirflowPlugin):
    name = "etna_plugin"

    hooks = [GitHook]
    appbuilder_menu_items = [appbuilder_mitem_toplevel]

    def on_load(*args, **kwargs):
        # ... perform Plugin boot actions
        logging.error("here I have been loaded")