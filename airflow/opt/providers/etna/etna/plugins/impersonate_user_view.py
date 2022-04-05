from flask import flash

from wtforms import SelectField
from wtforms.validators import InputRequired
from flask_appbuilder.fieldwidgets import Select2Widget
from flask_appbuilder.forms import DynamicForm


class ImpersonateUserForm(DynamicForm):
    username = SelectField(label='Username',
        validate_choice=False, # Seems to be required to work with Airflow, since the choices are injected
        description=('User whose permissions to impersonate'),
        validators = [InputRequired()],
        widget=Select2Widget())


try:
    from flask_appbuilder import SimpleFormView
    from flask_appbuilder.security.decorators import has_access

    class ImpersonateUserView(SimpleFormView):
        form = ImpersonateUserForm
        form_title = 'User whose permissions you want to impersonate'

        @has_access
        def form_get(self, form):
            def sort_user_email(u):
                return u.email

            users = self.appbuilder.sm.get_all_users()
            users.sort(key=sort_user_email)
            form.username.choices = [(u.email, u.email) for u in users]

        @has_access
        def form_post(self, form):
            # post process form
            subject = self.appbuilder.sm.find_user(email=form.username.data)

            if subject is None:
                raise ValueError(f"{form.username.data} not found.")

            subject_roles = self.appbuilder.sm.get_user_roles(subject)

            impersonator = self.appbuilder.sm.current_user

            impersonator.roles = subject_roles
            self.appbuilder.sm.update_user(impersonator)

            flash(f"Done -- you have the same Airflow permissions as {form.username.data}. Logout to reset your permissions.", 'info')

    impersonate_user_view = ImpersonateUserView()

except (ImportError, ModuleNotFoundError) as e:
    impersonate_user_view = None