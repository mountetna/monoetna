import contextlib
import functools
import logging
import os
from random import Random
from unittest.mock import patch

import jinja2
import pytest
from airflow import settings
from airflow.jobs.base_job import BaseJob
from airflow.models import (
    DagRun,
    TaskInstance,
    DagTag,
    DagModel,
    SlaMiss,
    Pool,
    RenderedTaskInstanceFields,
    errors,
    XCom,
    Log,
    TaskFail,
    TaskReschedule,
    Connection,
)
from airflow.models.dagbag import DagBag
from airflow.models.dagcode import DagCode
from airflow.models.serialized_dag import SerializedDagModel
from airflow.models.variable import Variable
from airflow.utils.db import add_default_pool_if_not_exists, create_default_connections
from airflow.utils.session import create_session
from airflow.www.app import purge_cached_app, create_app
from flask_appbuilder.security.sqla.models import User, PermissionView
from requests import Session
from sqlalchemy.orm import scoped_session, sessionmaker
from urllib3 import connectionpool
from vcr.patch import CassettePatcherBuilder
from vcr.stubs import VCRFakeSocket

from etna.hooks.etna import EtnaHook


class NotSoRandom(Random):
    def __init__(self, *args):
        super(NotSoRandom, self).__init__(0)


@pytest.fixture()
def session_a():
    from airflow.settings import engine

    Session = scoped_session(
        sessionmaker(
            autocommit=False,
            autoflush=False,
            bind=engine,
            expire_on_commit=False,
        )
    )
    session = Session()

    try:
        yield session
    finally:
        session.close()


@pytest.fixture()
def session_b():
    from airflow.settings import engine

    Session = scoped_session(
        sessionmaker(
            autocommit=False,
            autoflush=False,
            bind=engine,
            expire_on_commit=False,
        )
    )
    session = Session()

    try:
        yield session
    finally:
        session.close()


@pytest.fixture(scope="module")
def vcr_config():
    vcr_log = logging.getLogger("vcr")
    # change this to logging level to increase debugging, but it is very very verbose by default
    # so we keep it at warning.
    vcr_log.setLevel(logging.WARNING)

    # This greatly slows down connections in our tests, but unfortunately vcrpy does not
    # support keep alive correctly, meaning that tests that happen to last longer than the
    # keep alive configured for a server will fail with this error.
    class ForceNeverReuseConnectionSession(Session):
        def __init__(self, create_session, auth):
            self.create_session = create_session
            super(ForceNeverReuseConnectionSession, self).__init__()
            self.auth = auth

        def __call__(self):
            return self

        def request(self, *args, **kwds):
            with self.create_session() as session:
                return session.request(*args, **kwds)

    get_client = EtnaHook.get_client

    @contextlib.contextmanager
    def _get_client(self, auth):
        yield ForceNeverReuseConnectionSession(lambda: get_client(self, auth), auth)

    EtnaHook.get_client = _get_client

    # Authorization in any vcr is safely masked
    return {"filter_headers": [("authorization", "XXX-Auth")]}


# def match_body(r1: BaseRequest, r2: BaseRequest):
#     return r1.re
#
# def pytest_recording_configure(config, vcr):
#     vcr.register_matcher("body", match_body)


@pytest.fixture()
def rsa_etna_connection(session):
    session.query(Connection).filter(
        Connection.conn_id == "rsa_etna_connection"
    ).delete()

    conn = Connection(
        conn_id="rsa_etna_connection",
        conn_type="etna",
        login=os.environ.get("AIRFLOW_ETNA_TEST_EMAIL", "zachary.collins@ucsf.edu"),
        # This RSA key was cycled out after recording the test.  To record e2e again, try providing a new key, and replace it after recording.
        password=os.environ.get(
            "AIRFLOW_ETNA_TEST_RSA",
            "-----BEGIN PRIVATE KEY----- MIIEvgIBADANBgkqhkiG9w0BAQEFAASCBKgwggSkAgEAAoIBAQCX5/oljASRyfxf OG6uST61p5N36U7vsRgs0c2ZwsCfDNhl7drxglIVJRMshg5ZJazxmpMp15yb6u1b Y8iUy1ebPFn1BquuJIvNn2eRhV/zCbfETlbJXGT700jXemrdC8XFe9CCEg/J5lHr PkIEKKxY57LcpmHBHzo4f27Sn1PJzmdYlekSu1gLLbKxUmo2IexKpXmwTjTPViaa JFUsLmC07uzWrXcDo3JTjzz3Pk3KM8c+s9XLAl+ovrqO28vPghB5KzlSW8BCgQDc m0EdljkQIf3WdSNrhpIVUtWYLk29CJExgvY9Ei+8Ulw7gksN7+oDxcQefMNDgxQQ 4WPmClALAgMBAAECggEABy0WTB/JN3nrSjRIRkN/iuVXuhpzeC9NjRB8Pf9NSjY5 IteRuEcHyafut/O9ScjV2rQKr7dX1qXKgL6+Awl4IgU/2qtuANQJJrWZFu7OEZUr 8UIiJ3EN9DePAV7vHXIo7aNjvkFMLaWLySkvxTKGscyATpwtkgn/nhunCJwuQSJE baVCUEXN9OlMOcsTcuVmv4wEeLN+TRfMeMaLyv+WNDI29BZzuzny3qvr52ooYjjr 9HdYF9L+LM+LkuolMQIIy7jcAV0Y1RZ3F+HBQqc7qFOa+zFEQC9byr0+EMryoDU6 ZurPErgr2d3Wqo4g920tIn8TOannU6YXXNVIBrCFAQKBgQDDq1TwuvHFLuE/vXsP PGZUFMehcQS+f2Nm1iP8UIJbhmCCb0gxx2LlTxEiHKSq3AyA2AK3A5VOFQLmZUjM J49Fcv9pkPz2whYqCBk3qqwWCqBkYFmzQmiSaYyE19u7s/FTm67OSrguwPO3Euux /g0UVP4EXzlT0V9dAERYFJxu+wKBgQDGvlAwHOd8zngDkKyWBD8YyyL4RT5MluK+ ZShSHYrYjHRbVG05MIl5Hw1v6VA6DSLLlA//weD/EkS6KRnvhPP3h7AiQkV39OYx l5R25KbNuX6uETVIWcUHb34lQTsipI4YMtWlLKA9q9mjzZEDw/egIN8XAP08DSTK Aa0ujJaWMQKBgQClbwyH5GdZogNMEvYisZyK5m7KrnWmYqo2XkNapu8wVvLuFQxj GgMhgbIotzL6SsY/gWL6PYtU0yr6hRQBmEjoHQyZwr4+G2cF7obzq9eHY0Cs3VG5 4CHt+FOYVbEwiDk3yV8Ih+All3n3hYXFndiNIjcKl0Au/8yzIvClz/dbVQKBgQCP f+q2UqhyXUIakOOMjhRg+ouNZ7HL60Zc4v1yDRKruP5q01Lp8DnS0rEJFRVwVPvC sm265WpnwfEN2Y94ei8Nk1OB6QfvzUxIkoIINqCZ+k2VsacfTnINJFuY2riwEtDm eA367XXmEadbtpn2dhDd9d4e5f/y1Cq0EPHSooA4gQKBgHAx9MMGEIFc2z3Hc9AW E6S1+duQLnv5LBAxK4xgXwNchFNwa7QabXWo3nArKC5P/ORshMZZly5zLXRCX5Cy cf7CTPJJFxVrS7arUBPL4YzqsDzkzmSrpURc+fvZEaiWTnBWhYJPmZjwWl3sKqyY 5hIZhQjPLdvWuYJv4vE5qcMk -----END PRIVATE KEY----- ",
        ),
        schema="rsa",
        host=".ucsf.edu",
    )

    session.add(conn)
    session.commit()
    return conn


@pytest.fixture()
def token_etna_connection(session):
    session.query(Connection).filter(
        Connection.conn_id == "rsa_etna_connection"
    ).delete()

    conn = Connection(
        conn_id="rsa_etna_connection",
        conn_type="etna",
        login=os.environ.get("AIRFLOW_ETNA_TEST_EMAIL", ""),
        # Unlike the rsa, there is no client side logic that depends on this, it is safe to leave it empty
        # outside of recording new tests with new tokens.  It will be filtered from recordings.
        password=os.environ.get("AIRFLOW_ETNA_TEST_TOKEN", ""),
        schema="token",
        host="",
    )

    session.add(conn)
    session.commit()
    return conn


@pytest.fixture()
def https_git_connection(session):
    session.query(Connection).filter(
        Connection.conn_id == "https_git_connection"
    ).delete()

    conn = Connection(
        conn_id="https_git_connection",
        conn_type="git",
        login=os.environ.get("AIRFLOW_GIT_TEST_USER", ""),
        password=os.environ.get("AIRFLOW_GIT_TEST_PASSWORD", ""),
        schema="https",
    )

    session.add(conn)
    session.commit()
    return conn


@pytest.fixture()
def ssh_git_connection(session):
    password = os.environ.get("AIRFLOW_GIT_TEST_PK", "")
    session.query(Connection).filter(
        Connection.conn_id == "ssh_git_connection"
    ).delete()

    conn = Connection(
        conn_id="ssh_git_connection",
        conn_type="git",
        login=os.environ.get("AIRFLOW_GIT_TEST_USER", ""),
        password=password,
        schema="ssh",
    )

    session.add(conn)
    session.commit()
    return conn


@pytest.fixture()
def reset_environment():
    """
    Resets env variables.
    """
    init_env = os.environ.copy()
    yield
    changed_env = os.environ
    for key in changed_env:
        if key not in init_env:
            del os.environ[key]
        else:
            os.environ[key] = init_env[key]


@pytest.fixture(scope="module")
def session():
    settings.configure_orm()
    yield settings.Session


@pytest.fixture()
def reset_db():
    from airflow.utils import db

    db.resetdb()
    yield


@pytest.fixture(scope="module")
def examples_dag_bag(session):
    DagBag(include_examples=True).sync_to_db()
    dag_bag = DagBag(include_examples=True, read_dags_from_db=True)
    session.commit()
    yield dag_bag


def dont_initialize_flask_app_submodules(_func=None, *, skip_all_except=None):
    if not skip_all_except:
        skip_all_except = []

    def decorator_dont_initialize_flask_app_submodules(f):
        def no_op(*args, **kwargs):
            pass

        methods = [
            "init_api_experimental_auth",
            "init_flash_views",
            "init_appbuilder_links",
            "init_appbuilder_views",
            "init_plugins",
            "init_connection_form",
            "init_error_handlers",
            "init_api_connexion",
            "init_api_experimental",
            "sync_appbuilder_roles",
            "init_jinja_globals",
            "init_xframe_protection",
            "init_appbuilder",
        ]

        @functools.wraps(f)
        def func(*args, **kwargs):

            for method in methods:
                if method not in skip_all_except:
                    patcher = patch(f"airflow.www.app.{method}", no_op)
                    patcher.start()
            purge_cached_app()
            result = f(*args, **kwargs)
            patch.stopall()
            purge_cached_app()

            return result

        return func

    if _func is None:
        return decorator_dont_initialize_flask_app_submodules
    else:
        return decorator_dont_initialize_flask_app_submodules(_func)


@pytest.fixture(scope="module")
def app(examples_dag_bag):
    @dont_initialize_flask_app_submodules(
        skip_all_except=[
            "init_api_connexion",
            "init_appbuilder",
            "init_appbuilder_links",
            "init_appbuilder_views",
            "init_flash_views",
            "init_jinja_globals",
            "init_plugins",
        ]
    )
    def factory():
        return create_app(testing=True)

    app = factory()
    app.config["WTF_CSRF_ENABLED"] = False
    app.dag_bag = examples_dag_bag
    app.jinja_env.undefined = jinja2.StrictUndefined

    security_manager = app.appbuilder.sm

    test_users = [
        {
            "username": "test_admin",
            "first_name": "test_admin_first_name",
            "last_name": "test_admin_last_name",
            "email": "test_admin@fab.org",
            "role": security_manager.find_role("Admin"),
            "password": "test_admin_password",
        },
        {
            "username": "test_user",
            "first_name": "test_user_first_name",
            "last_name": "test_user_last_name",
            "email": "test_user@fab.org",
            "role": security_manager.find_role("User"),
            "password": "test_user_password",
        },
        {
            "username": "test_viewer",
            "first_name": "test_viewer_first_name",
            "last_name": "test_viewer_last_name",
            "email": "test_viewer@fab.org",
            "role": security_manager.find_role("Viewer"),
            "password": "test_viewer_password",
        },
        {
            "username": "test_project_viewer",
            "first_name": "test_project_viewer_first_name",
            "last_name": "test_project_viewer_last_name",
            "email": "test_project_viewer@fab.org",
            "role": security_manager.sync_project_role("test"),
            "password": "test_project_viewer_password",
        },
    ]

    for user_dict in test_users:
        if not security_manager.find_user(username=user_dict["username"]):
            security_manager.add_user(**user_dict)

    yield app


@pytest.fixture()
def admin_client(app):
    return app.test_client()


@pytest.fixture()
def reset_admin_role(app):
    yield
    user = app.appbuilder.sm.find_user(username='test_admin')
    admin_role = app.appbuilder.sm.find_role('Admin')
    user.roles = [admin_role]
    app.appbuilder.sm.update_user(user)


def clear_users():
    with create_session() as session:
        session.query(User).delete()
        session.query(PermissionView).delete()


def clear_db_runs():
    with create_session() as session:
        # session.query(TriggererJob).delete()
        # session.query(Trigger).delete()
        session.query(DagRun).delete()
        session.query(TaskInstance).delete()


def clear_db_dags():
    with create_session() as session:
        session.query(DagTag).delete()
        session.query(DagModel).delete()


def clear_db_serialized_dags():
    with create_session() as session:
        session.query(SerializedDagModel).delete()


def clear_db_sla_miss():
    with create_session() as session:
        session.query(SlaMiss).delete()


def clear_db_pools():
    with create_session() as session:
        session.query(Pool).delete()
        add_default_pool_if_not_exists(session)


def clear_db_connections(add_default_connections_back=True):
    with create_session() as session:
        session.query(Connection).delete()
        if add_default_connections_back:
            create_default_connections(session)


def clear_db_variables():
    with create_session() as session:
        session.query(Variable).delete()


def clear_db_dag_code():
    with create_session() as session:
        session.query(DagCode).delete()


def set_default_pool_slots(slots):
    with create_session() as session:
        default_pool = Pool.get_default_pool(session)
        default_pool.slots = slots


def clear_rendered_ti_fields():
    with create_session() as session:
        session.query(RenderedTaskInstanceFields).delete()


def clear_db_import_errors():
    with create_session() as session:
        session.query(errors.ImportError).delete()


def clear_db_xcom():
    with create_session() as session:
        session.query(XCom).delete()


def clear_db_logs():
    with create_session() as session:
        session.query(Log).delete()


def clear_db_jobs():
    with create_session() as session:
        session.query(BaseJob).delete()


def clear_db_task_fail():
    with create_session() as session:
        session.query(TaskFail).delete()


def clear_db_task_reschedule():
    with create_session() as session:
        session.query(TaskReschedule).delete()
