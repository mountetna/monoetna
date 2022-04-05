from airflow.plugins_manager import AirflowPlugin
from flask import Blueprint

import logging


ETNA_DOCS_MENU_TITLE = "Etna Airflow Docs"
ETNA_DOCS_STATIC_PATH = "/static/docs"

# Creating flask appbuilder Menu Items
appbuilder_mitem_toplevel = {
    "name": ETNA_DOCS_MENU_TITLE,
    "href": f"{ETNA_DOCS_STATIC_PATH}/index.html",
}

# Creating a flask blueprint to integrate the doc static folder
docs_blueprint = Blueprint(
    "etna_plugin_docs_blueprint",
    __name__,
    static_folder='/opt/airflow/providers/etna/docs',
    static_url_path=ETNA_DOCS_STATIC_PATH)


class EtnaDocsPlugin(AirflowPlugin):
    name = "etna_docs_plugin"

    appbuilder_menu_items = [appbuilder_mitem_toplevel]
    flask_blueprints = [docs_blueprint]
