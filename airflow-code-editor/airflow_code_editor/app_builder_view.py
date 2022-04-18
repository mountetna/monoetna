#!/usr/bin/env python
#
#   Copyright 2019 Andrea Bonomi <andrea.bonomi@gmail.com>
#   Modifications included by DSCOLabs engineering team
#   University of San Francisco, CA
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the Licens
#

from flask import request
from flask_appbuilder import BaseView, expose
from airflow_code_editor.code_editor_view import AbstractCodeEditorView
from airflow_code_editor.commons import (
    ROUTE,
    MENU_CATEGORY,
    MENU_LABEL,
    JS_FILES,
    VERSION,
)

__all__ = ["appbuilder_view"]

try:
    from airflow.www import auth
    from airflow.security import permissions

    PERMISSIONS = [
        (permissions.ACTION_CAN_READ, permissions.RESOURCE_WEBSITE),
    ]

    class AppBuilderCodeEditorView(BaseView, AbstractCodeEditorView):
        route_base = ROUTE
        base_permissions = ["can_list", "can_create", "menu_acccess"]

        @expose("/")
        @auth.has_access(PERMISSIONS)
        def list(self):
            return self._index()

        @expose("/repo", methods=["POST"])
        @auth.has_access(PERMISSIONS)
        def repo_base(self, path=None):
            return self._git_repo(path)

        @expose("/repo/<path:path>", methods=["GET", "HEAD", "POST"])
        @auth.has_access(PERMISSIONS)
        def repo(self, path=None):
            return self._git_repo(path)

        @expose("/files/<path:path>", methods=["POST"])
        @auth.has_access(PERMISSIONS)
        def save(self, path=None):
            return self._save(path)

        @expose("/uploads/<path:path>", methods=["POST"])
        @auth.has_access(PERMISSIONS)
        def upload(self, path=None):
            return self._save(path, upload=True)

        @expose("/files/<path:path>", methods=["GET"])
        @auth.has_access(PERMISSIONS)
        def load(self, path=None):
            return self._load(path)

        @expose("/format", methods=["POST"])
        @auth.has_access(PERMISSIONS)
        def format(self):
            return self._format()

        @expose("/tree", methods=["GET"])
        @auth.has_access(PERMISSIONS)
        def tree_base(self, path=None):
            return self._tree(path, args=request.args)

        @expose("/tree/<path:path>", methods=["GET"])
        @auth.has_access(PERMISSIONS)
        def tree(self, path=None):
            return self._tree(path, args=request.args)

        @expose("/ping", methods=["GET"])
        @auth.has_access(PERMISSIONS)
        def ping(self):
            return self._ping()

        def _render(self, template, *args, **kargs):
            return self.render_template(
                template + "_appbuilder.html",
                airflow_major_version=self.airflow_major_version,
                js_files=JS_FILES,
                version=VERSION,
                *args,
                **kargs
            )

except (ImportError, ModuleNotFoundError):
    # We do not support older airflow versions.
    raise


appbuilder_code_editor_view = AppBuilderCodeEditorView()
appbuilder_view = {
    "category": MENU_CATEGORY,
    "name": MENU_LABEL,
    "view": appbuilder_code_editor_view,
}
