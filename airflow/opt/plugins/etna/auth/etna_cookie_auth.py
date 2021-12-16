import logging
from typing import List
from urllib.parse import urlencode

import jwt
from airflow.models.dag import DagModel, DagTag
from airflow.security import permissions
from airflow.security.permissions import resource_name_for_dag
from airflow.www.security import AirflowSecurityManager
from airflow.www.views import CustomUserRemoteUserModelView
from flask import redirect, request
from flask_appbuilder import expose
from flask_appbuilder.security.sqla.models import Role, User
from flask_appbuilder.security.views import AuthView
from flask_login import login_user
from werkzeug.wrappers import Response as WerkzeugResponse

from etna.auth.etna_user import EtnaUser


def deserialize_etna_user(payload: str, algo: str, key: str) -> EtnaUser:
    user_json = jwt.decode(
        payload,
        key=key,
        verify=True,
        algorithms=[algo],
    )

    return EtnaUser.from_dict(user_json)


class EtnaCookieAuthView(AuthView):
    login_template = ""

    @expose("/login/")
    def login(self) -> WerkzeugResponse:
        cookie_name = self.appbuilder.app.config["ETNA_AUTH_COOKIE_NAME"]
        cookie = request.cookies.get(cookie_name)
        user = None

        logging.info(str(cookie))

        if cookie:
            user = self.apply_etna_user(cookie, self.appbuilder.sm)

        if not user:
            return redirect(
                self.appbuilder.app.config["ETNA_AUTH_JANUS_URL"]
                + "/login?"
                + urlencode(dict(referer=request.url))
            )

        return redirect(self.appbuilder.get_url_for_index)

    def apply_etna_user(self, cookie: str, security_manager: "EtnaSecurityManager"):
        token_algo = self.appbuilder.app.config["ETNA_AUTH_TOKEN_ALGO"]
        public_key = self.appbuilder.app.config["ETNA_AUTH_PUBLIC_KEY"]
        etna_user = deserialize_etna_user(cookie, token_algo, public_key)

        user = security_manager.find_user(username=etna_user.username)
        if user is None:
            user = security_manager.add_user(
                username=etna_user.username,
                first_name=etna_user.name,
                last_name="-",
                email=etna_user.email,
                role=[],
            )

        security_manager.sync_user(user, etna_user)
        login_user(user)
        security_manager.update_user_auth_stat(user)
        return user


def role_name_by_project(project_name: str, role_name: str) -> str:
    return f"{project_name}_{role_name}"


class EtnaSecurityManager(AirflowSecurityManager):
    authremoteuserview = EtnaCookieAuthView

    def sync_user(self, user: User, etna_user: EtnaUser) -> User:
        user.roles = self.roles_for_etna_user(etna_user)
        self.update_user(user)
        return user

    def roles_for_etna_user(self, etna_user: EtnaUser) -> List[Role]:
        roles: List[Role] = []

        if etna_user.is_superuser:
            roles.append(self.find_role("Admin"))
        elif etna_user.is_supereditor:
            roles.append(self.find_role("User"))
        elif etna_user.is_superviewer:
            roles.append(self.find_role("Viewer"))

        for project in etna_user.projects:
            roles.append(self.sync_project_role(project))

        return roles

    def sync_project_role(self, project_name: str) -> Role:
        role_name = role_name_by_project(project_name, "viewer")
        role = self.find_role(role_name)

        if not role:
            role = self.add_role(role_name)

        dags_query = self.get_session.query(DagModel)
        cond = DagModel.tags.any(DagTag.name == f"project_{project_name}")
        dags_query = dags_query.filter(cond)

        for dag in dags_query.all():
            perm_view = self.add_permission_view_menu(
                permissions.ACTION_CAN_READ,
                resource_name_for_dag(dag.dag_id),
            )
            self.add_permission_role(role, perm_view)

        return role
