import re
from datetime import datetime

from airflow import DAG
from airflow.executors.debug_executor import DebugExecutor

from etna.hooks.etna import Folder, File
from etna.etls.metis import filter_by_record_directory


def run_dag(dag: DAG, execution_date: datetime, end_date: datetime):
    dag.run(
        executor=DebugExecutor(),
        start_date=execution_date,
        end_date=end_date,
        verbose=True,
    )

def test_filter_by_root_directory():
    result = filter_by_record_directory([
        Folder(folder_path="bulk_RNASeq/raw/abcdef"),
        File(file_path="bulk_RNASeq/raw/abcdef/myfile.txt"),
        File(file_path="bulk_RNASeq/raw/abcdef2/myfile.txt"),
    ], re.compile(r'bulk_RNASeq/raw/[^/]*'), 'model')

    assert result[0].root_path == 'bulk_RNASeq/raw/abcdef'
    assert result[0].record_name == 'abcdef'
    assert result[0].match_subpath == ''
    assert result[0].match_full_path == 'bulk_RNASeq/raw/abcdef'
    assert result[0].folder_path == 'bulk_RNASeq/raw/abcdef'
    assert result[0].match_file is None
    assert result[0].match_folder is not None

    assert result[1].root_path == 'bulk_RNASeq/raw/abcdef2'
    assert result[1].record_name == 'abcdef2'
    assert result[1].match_subpath == 'myfile.txt'
    assert result[1].match_full_path == 'bulk_RNASeq/raw/abcdef2/myfile.txt'
    assert result[1].folder_path == 'bulk_RNASeq/raw/abcdef2'
    assert result[1].match_file is not None
    assert result[1].match_folder is None


#
# @mock.patch("tempfile._Random", NotSoRandom)
# def test_metis_files_etl_e2e(token_etna_connection: Connection):
#     hook = EtnaHook(token_etna_connection.conn_id)
#
#     start = datetime(2021, 1, 1)
#
#     def t(i=0, days=0):
#         return (start + timedelta(minutes=i, days=days)).replace(tzinfo=utc)
#
#     @metis_etl('mvir1', 'data', datetime(2020, 1, 1), 1, hook=hook)
#     def load_metis_files_etl_dag(tail_folders, tail_files):
#         @task
#         def find_rna_seq_record_folders(folders: List[Folder], files: List[File]):
#             matching_record_folders = filter_by_record_directory(folders + files, re.compile(r'^bulk_RNASeq/raw/[^/]*'))
#             with hook.magma() as magma:
#                 matching_record_folders = filter_by_exists_in_timur(magma, matching_record_folders, 'rna_seq')
#             with hook.metis() as metis:
#                 return pickled(list_contents_of_matches(metis, matching_record_folders))
#
#         @link('rna_seq', 'raw_fastq_files', dry_run=True, hook=hook)
#         def link_raw_fastq_files(matches: List[Tuple[MatchedAtRoot, List[File]]]):
#             for match, files in matches:
#                 yield match, [f for f in files if re.compile(r'.*\.fastq\.gz$').match(f.file_path)]
#
#         rna_seq_record_folders = find_rna_seq_record_folders(tail_folders, tail_files)
#         link_raw_fastq_files(rna_seq_record_folders)
#
#         cat_in_container = run_in_container('cat_in_container', 'polyphemus_app', ['cat', '/inputs/batch'],
#                                             output_json=List[File], docker_base_url="http://localhost:8085")
#         cat_in_container['/inputs/batch'] = tail_folders
#
#     load_metis_files_etl_dag: DAG
#
#     run_dag(load_metis_files_etl_dag, t(), t(4))
#     folders = XCom.get_one(dag_id=load_metis_files_etl_dag.dag_id, task_id='tail_folders', execution_date=t())
#     assert len(folders) > 0
#
#     run_dag(load_metis_files_etl_dag, t(6), t(9))
#     folders = XCom.get_one(dag_id=load_metis_files_etl_dag.dag_id, task_id='tail_folders', execution_date=t(6))
#     assert len(folders) == 0
