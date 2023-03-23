import re
import pytz
from datetime import timedelta, datetime
from dateutil import parser
from unittest import mock
from pytest import raises

from airflow import DAG
from airflow.decorators import task
from airflow.models.xcom_arg import XComArg
from airflow.models import Variable, Connection
from airflow.utils.state import State


from providers.etna.etna.etls.universal_metis_linker import UniversalMetisLinker
from providers.etna.etna.etls.box import BoxEtlHelpers
from etna import system_dag

from .test_metis_files_etl import run_dag, get_all_results

mock_etna_hook = mock.Mock()
mock_polyphemus = mock.Mock()
mock_polyphemus.return_value.__enter__ = mock_polyphemus
mock_polyphemus.return_value.__exit__ = mock_polyphemus
mock_etna_hook.polyphemus = mock_polyphemus
from etna.hooks.etna import EtnaHook

@mock.patch.object(EtnaHook, 'for_project', return_value=mock_etna_hook)
def test_universal_metis_linker(token_etna_connection: Connection):

    @system_dag(timedelta(minutes=2))
    def test_univ_linker():
        @task
        def run_linker():
            #linker = UniversalMetisLinker()
            #linker.hook.etna_conn_id = token_etna_connection.conn_id

            pass

            #linker.run()

        run_linker()

    start_date = parser.parse("2022-01-01 00:00:00 +0000")
    end_date = start_date + timedelta(days=1, minutes=1)
    test_univ_linker.clear()
    run_dag(test_univ_linker, start_date, end_date)
