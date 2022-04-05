from flask import redirect, flash
import airflow
from functools import wraps

import logging

from wtforms import SelectField
from wtforms.validators import InputRequired
from flask_appbuilder.fieldwidgets import Select2Widget
from flask_appbuilder.forms import DynamicForm


class ImpersonateUserForm(DynamicForm):
    username = SelectField(label='Username',
        coerce=int,
        validate_choice=False, # Seems to be required to work with Airflow, since the choices are injected
        description=('User whose permissions to impersonate'),
        validators = [InputRequired()],
        widget=Select2Widget())


try:
    from flask_appbuilder import SimpleFormView, expose
    from flask_appbuilder.security.decorators import has_access

    def superuser_required(func):
        # when airflow loads plugins, login is still None.
        @wraps(func)
        def func_wrapper(*args, **kwargs):
            if airflow.login and airflow.login.current_user.is_superuser():
                return func(*args, **kwargs)
            return redirect(args[0].appbuilder.get_url_for_index)

        return func_wrapper

    class ImpersonateUserView(SimpleFormView):
        form = ImpersonateUserForm
        form_title = 'User whose permissions you want to impersonate'
        message = 'Submitted'

        @has_access
        def form_get(self, form):
            users = self.appbuilder.sm.get_all_users()
            form.username.choices = [(u.id, u.email) for u in users]

        @has_access
        def form_post(self, form):
            # post process form
            flash(self.message, 'info')
        # default_view = "index"

        # @expose("/")
        # @superuser_required
        # def index(self):
        #     return self.render(
        #         "impersonate_user.html",
        #         airflow_major_version=self.airflow_major_version,
        #         js_files=["impersonate_helpers.js"],
        #         version="0.0.1",
        #     )

        # @expose("/submit", methods=["POST"])
        # @superuser_required
        # def impersonate(self, path=None):
        #     # Do something here to switch the user's permissions
        #     return redirect(self.appbuilder.get_url_for_index)

    impersonate_user_view = ImpersonateUserView()

except (ImportError, ModuleNotFoundError) as e:
    logging.error("error!")
    logging.error(e)
    impersonate_user_view = None