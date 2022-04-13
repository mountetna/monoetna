import re
import os
import os.path

from airflow import configuration
from flask import flash

from wtforms import StringField
from wtforms.validators import Regexp, InputRequired
from flask_appbuilder.forms import DynamicForm
from flask_appbuilder import SimpleFormView
from flask_appbuilder.security.decorators import has_access


class AddProjectForm(DynamicForm):
    project_name = StringField(
        label='Project name',
        description="Short name of the project to be added",
        validators=[InputRequired(), Regexp(re.compile("(?!pg_)[a-z][a-z0-9]*(_[a-z][a-z0-9]*)*"))]
    )

class AddProjectView(SimpleFormView):
    form = AddProjectForm
    form_title = 'Project you want to add to the dag code editor'

    @has_access
    def form_get(self, form):
        pass

    @has_access
    def form_post(self, form):
        dags_folder = configuration.conf.get('core', 'dags_folder')
        os.mkdir(os.path.join(dags_folder, form.project_name.data))
        flash(f"Done -- a project dags directory has been created for {form.project_name.data}", 'info')

add_project_view = AddProjectView()
