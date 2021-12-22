import functools
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

class NotSoRandom(Random):
    def __init__(self, *args):
        super(NotSoRandom, self).__init__(0)

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
            "init_permanent_session",
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
    ]

    for user_dict in test_users:
        if not security_manager.find_user(username=user_dict["username"]):
            security_manager.add_user(**user_dict)

    yield app


@pytest.fixture()
def admin_client(app):
    return app.test_client()


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
