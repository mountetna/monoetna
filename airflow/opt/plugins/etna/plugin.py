from airflow.plugins_manager import AirflowPlugin
from flask import Blueprint

from etna.hooks.git import GitHook

import logging


ETNA_DOCS_MENU_TITLE = "Etna Airflow Docs"

# Creating flask appbuilder Menu Items
appbuilder_mitem_toplevel = {
    "name": ETNA_DOCS_MENU_TITLE,
    "href": "/etnadocs",
}

# Creating a flask blueprint to integrate the doc static folder
bp = Blueprint(
    "etna_docs_blueprint", __name__,
    static_folder='/opt/airflow/providers/etna/docs',
    static_url_path='/etnadocs')

class EtnaPlugin(AirflowPlugin):
    name = "etna_plugin"

    hooks = [GitHook]
    appbuilder_menu_items = [appbuilder_mitem_toplevel]
    flask_blueprints = [bp]

    def on_load(*args, **kwargs):
        # ... perform Plugin boot actions
        logging.error("here I have been loaded")