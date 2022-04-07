from airflow.security import permissions

from etna.plugins.etna_docs_plugin import ETNA_DOCS_MENU_TITLE
from airflow_code_editor.commons import MENU_LABEL as CODE_EDITOR_MENU_LABEL


def test_user_project_permissions(app):
    def to_text(action: str, element: str):
        return f"{action.replace('_', ' ')} on {element}"

    sm = app.appbuilder.sm
    role = sm.sync_project_role('test')
    perm_strs = [str(p) for p in role.permissions]

    assert(to_text(permissions.ACTION_CAN_ACCESS_MENU, ETNA_DOCS_MENU_TITLE) in perm_strs)
    assert(to_text(permissions.ACTION_CAN_ACCESS_MENU, CODE_EDITOR_MENU_LABEL) in perm_strs)
    assert(to_text(permissions.ACTION_CAN_READ, permissions.RESOURCE_CONNECTION) not in perm_strs)
    assert(to_text(permissions.ACTION_CAN_READ, permissions.RESOURCE_DAG) not in perm_strs) # Cannot read all DAGs