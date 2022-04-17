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
import importlib
import importlib.machinery
import importlib.util
import os
import os.path
import logging
import mimetypes

import sys
import tempfile

from airflow.configuration import conf
from airflow.models.dagbag import DagBag
from airflow.utils.timeout import timeout
from flask import abort, request, send_file
from flask_wtf.csrf import generate_csrf
from airflow.version import version
from airflow_code_editor.commons import HTTP_404_NOT_FOUND
from airflow_code_editor.tree import get_tree
from airflow_code_editor.utils import (
    get_plugin_boolean_config,
    get_plugin_int_config,
    git_absolute_path,
    execute_git_command,
    error_message,
    normalize_path,
    prepare_api_response,
)


__all__ = ["AbstractCodeEditorView"]

AIRFLOW_MAJOR_VERSION = int(version.split(".")[0])

def clean_data(data):
    return data.replace("\r", "").rstrip()

class AbstractCodeEditorView(object):
    airflow_major_version = AIRFLOW_MAJOR_VERSION

    def _index(self):
        return self._render("index")

    def _save(self, path=None, upload=False):
        try:
            fullpath = git_absolute_path(path)
            if not upload:
                payload = request.get_json()
                data = payload['to']

                # Newline fix (remove cr)
                data = clean_data(data)
                os.makedirs(os.path.dirname(fullpath), exist_ok=True)

                # "Version" the code so as to prevent concurrent modifications.
                if os.path.exists(fullpath):
                    with open(fullpath, 'r') as f:
                        if clean_data(f.read()) != clean_data(payload['from']):
                            return prepare_api_response(path=normalize_path(path), error_message="The file has been modified since you last opened it!  Copy your changes and refresh the page.")

                if fullpath.endswith('.py'):
                    # Validate the file, write it as a temporary object and try importing it in the dag bag context.
                    tempdir = tempfile.mkdtemp()
                    with tempfile.NamedTemporaryFile(dir=tempdir, suffix='.py') as file:
                        file.write(data.encode('utf-8'))
                        file.write("\n")

                        try:
                            with timeout(conf.getfloat('core', 'DAGBAG_IMPORT_TIMEOUT'), error_message="Timeout importing module"):
                                mod_name = "airflow_code_editor_test_module"
                                loader = importlib.machinery.SourceFileLoader(mod_name, file.name)
                                spec = importlib.util.spec_from_loader(mod_name, loader)
                                new_module = importlib.util.module_from_spec(spec)
                                loader.exec_module(new_module)
                        except Exception as e:
                            return prepare_api_response(
                                path=normalize_path(path),
                                error_message=f"Problem loading module: {e}"
                            )

                with open(fullpath, "w") as f:
                    f.write(data)
                    f.write("\n")
            else:
                data = request.get_data()
                os.makedirs(os.path.dirname(fullpath), exist_ok=True)
                with open(fullpath, "wb") as f:
                    f.write(data)
            return prepare_api_response(path=normalize_path(path))
        except Exception as ex:
            logging.error(ex)
            return prepare_api_response(
                path=normalize_path(path),
                error_message="Error saving {path}: {message}".format(
                    path=path, message=error_message(ex)
                ),
            )

    def _git_repo(self, path):
        if request.method == "POST":
            return self._git_repo_post(path)
        else:
            return self._git_repo_get(path)

    def _git_repo_get(self, path):
        " Get a file from GIT (invoked by the HTTP GET method) "
        return execute_git_command(["cat-file", "-p", path])

    def _git_repo_post(self, path):
        " Execute a GIT command (invoked by the HTTP POST method) "
        git_args = request.json.get('args', [])
        return execute_git_command(git_args)

    def _load(self, path):
        " Send the contents of a file to the client "
        try:
            path = normalize_path(path)
            if path.startswith("~git/"):
                # Download git blob - path = '~git/<hash>/<name>'
                _, path, filename = path.split("/", 3)
                response = execute_git_command(["cat-file", "-p", path])
                response.headers["Content-Disposition"] = (
                    'attachment; filename="%s"' % filename
                )
                try:
                    content_type = mimetypes.guess_type(filename)[0]
                    if content_type:
                        response.headers["Content-Type"] = content_type
                except Exception:
                    pass
                return response
            else:
                # Download file
                fullpath = git_absolute_path(path)
                return send_file(fullpath, as_attachment=True)
        except Exception as ex:
            logging.error(ex)
            abort(HTTP_404_NOT_FOUND)

    def _format(self):
        " Format code "
        try:
            import black

            data = request.get_data(as_text=True)
            # Newline fix (remove cr)
            data = data.replace("\r", "").rstrip()
            mode = black.Mode(
                string_normalization=get_plugin_boolean_config("string_normalization"),
                line_length=get_plugin_int_config("line_length"),
            )
            data = black.format_str(src_contents=data, mode=mode)
            return prepare_api_response(data=data)
        except ImportError:
            return prepare_api_response(
                error_message="black dependency is not installed: to install black `pip install black`"
            )
        except Exception as ex:
            logging.error(ex)
            return prepare_api_response(
                error_message="Error formatting: {message}".format(
                    message=error_message(ex)
                )
            )

    def _tree(self, path, args = {}):
        return {'value': get_tree(path, args)}

    def _ping(self):
        return {'value': generate_csrf()}
