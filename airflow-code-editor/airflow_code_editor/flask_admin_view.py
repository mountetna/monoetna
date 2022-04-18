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
import airflow
from functools import wraps
from airflow_code_editor.code_editor_view import AbstractCodeEditorView
from airflow_code_editor.commons import (
    ROUTE,
    MENU_CATEGORY,
    MENU_LABEL,
    JS_FILES,
    VERSION,
)

__all__ = ["admin_view"]

# ############################################################################
# Flask Admin (Airflow < 2.0 and rbac = False)

# We do not support older airflow versions.
admin_view = None
