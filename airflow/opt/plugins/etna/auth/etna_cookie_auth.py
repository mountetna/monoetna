import os

import jwt
from flask_login import login_user, logout_user
from airflow.www.security import AirflowSecurityManager
from flask import flash, g, redirect, request
from flask_appbuilder import expose
from flask_appbuilder.security.views import AuthView
from werkzeug.wrappers import Response as WerkzeugResponse


class EtnaSecurityManager(AirflowSecurityManager):
  def auth_user_remote_user(self, username):
    pass

  pass


def deserialize_etna_user(payload: str, algo: str, key: str):
  try:
    user_json = jwt.decode(
      payload,
      key=key,
      verify=True,
      algorithms=[algo],
    )
  except jwt.DecodeError:
    pass


class EtnaCookieAuthView(AuthView):
  login_template = ""

  @expose("/login/")
  def login(self) -> WerkzeugResponse:
    cookie_name = self.appbuilder.app.config['ETNA_AUTH_COOKIE_NAME']
    cookie = request.cookies.get(cookie_name)
    if cookie:
      self._try_auth(cookie)

    return redirect(self.appbuilder.get_url_for_index)

  def _try_auth(self, cookie: str):
    token_algo = self.appbuilder.app.config['ETNA_AUTH_TOKEN_ALGO']
    public_key = self.appbuilder.app.config['ETNA_AUTH_PUBLIC_KEY']
    return deserialize_etna_user(cookie, token_algo, public_key)



    # username = request.environ.get("REMOTE_USER")
    # if g.user is not None and g.user.is_authenticated:
    #   return redirect(self.appbuilder.get_url_for_index)
    # if username:
    #   user = self.appbuilder.sm.auth_user_remote_user(username)
    #   if user is None:
    #     flash("invalid login", "warning")
    #   else:
    #     login_user(user)
    # else:
    #   flash("invalid login", "warning")
    # return redirect(self.appbuilder.get_url_for_index)
