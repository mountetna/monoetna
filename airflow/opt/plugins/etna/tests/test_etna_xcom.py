from datetime import datetime

from airflow.models import XCom
from airflow.utils.timezone import utc

from etna.hooks.etna import Folder
from etna.xcom.etna_xcom import pickled


def test_etna_xcom_pickled():
    now = datetime.now().replace(tzinfo=utc)
    XCom.set("key", pickled([Folder(folder_path="/abc")]), "task_id", "dag_id", now)
    assert XCom.get_one(execution_date=now, dag_id="dag_id", task_id="task_id") == [Folder(folder_path="/abc")]