import re
from io import StringIO
import pytz
from dataclasses import asdict
from datetime import timedelta, datetime
from dateutil import parser
from unittest import mock
from pytest import raises

from airflow import DAG
from airflow.decorators import task
from airflow.models.xcom_arg import XComArg
from airflow.models import Variable, Connection
from airflow.utils.state import State
from copy import deepcopy


from providers.etna.etna.etls.metis_linker import MetisLoaderConfig, MetisLinker, MetisLoaderError
from providers.etna.etna.etls.box import BoxEtlHelpers
from etna import system_dag
from mountetna import TailResultContainer, TailNode, Template, Attribute
from serde.json import from_json, to_json

def prep_tail(project_name, bucket_name, lines):
    container = TailResultContainer(bucket_name, project_name)
    for line in lines:
        container.add(TailNode(**line))

    return container.resolve_files(), container.resolve_folders(include_parents=True)


class TestMetisLinker:
    def test_matches_tails_to_file_scripts(self):
        # we have a tail
        tail = prep_tail('labors', 'pics', [
            {"type":"file","id":50,"parent_id":None,"node_name":"ignore.deceased.png","updated_at":"2023-08-03 22:39:17 +0000","file_hash":"0cc175b9c0f1b6a831c399e269772661","archive_id":None},
            {"type":"file","id":51,"parent_id":40,"node_name":"LABORS-LION-H2-C1.deceased.png","updated_at":"2023-11-11 22:39:17 +0000","file_hash":"8277e0910d750195b448797616e091ad","archive_id":None},
            {"type":"file","id":52,"parent_id":40,"node_name":"NO-ID.deceased.png","updated_at":"2023-11-11 22:39:17 +0000","file_hash":"8277e0910d750195b448797616e091ad","archive_id":None},
            {"type":"file","id":53,"parent_id":40,"node_name":"nonmatch_deceased.png","updated_at":"2023-11-11 22:39:17 +0000","file_hash":"8277e0910d750195b448797616e091ad","archive_id":None},
            {"type":"parent","id":40,"parent_id":None,"node_name":"victims","updated_at":"2023-08-03 22:39:17 +0000","file_hash":None,"archive_id":None},
            {"type":"parent","id":41,"parent_id":40,"node_name":"archived","updated_at":"2024-02-19 22:39:17 +0000","file_hash":None,"archive_id":None},
        ])

        # we have a config
        config = {
            "project_name": "labors",
            "config": {
                "bucket_name": "pics",
                "models": {
                    "victim": {
                        "scripts": [
                            {
                                "type": "file",
                                "folder_path": "victims",
                                "file_match": "*.deceased*.png",
                                "attribute_name": "photo_deceased"
                            }
                        ]
                    }
                }
            }
        }

        rules={
            'victim': r'^LABORS-LION-H\d+-C\d+$'
        }

        # given these we get an update
        update = MetisLoaderConfig(**config, rules=rules).update_for(tail)

        assert asdict(update)['revisions'] == {
            'victim': {
                'LABORS-LION-H2-C1': {
                    'photo_deceased': {
                        'original_filename': 'LABORS-LION-H2-C1.deceased.png',
                        'path': 'metis://labors/pics/victims/LABORS-LION-H2-C1.deceased.png'
                    }
                }
            }
        }

    def test_matches_tails_to_file_collection_scripts(self):
        # we have a tail
        tail = prep_tail('labors', 'pics', [
            {"type":"file","id":50,"parent_id":None,"node_name":"ignore.png","updated_at":"2023-08-03 22:39:17 +0000","file_hash":"0cc175b9c0f1b6a831c399e269772661","archive_id":None},
            {"type":"file","id":51,"parent_id":41,"node_name":"family_photo.1.png","updated_at":"2023-11-11 22:39:17 +0000","file_hash":"8277e0910d750195b448797616e091ad","archive_id":None},
            {"type":"file","id":52,"parent_id":41,"node_name":"family_photo.2.png","updated_at":"2023-11-11 22:39:17 +0000","file_hash":"9277e0910d750195b448797616e091ad","archive_id":None},
            {"type":"parent","id":40,"parent_id":None,"node_name":"family","updated_at":"2023-08-03 22:39:17 +0000","file_hash":None,"archive_id":None},
            {"type":"parent","id":41,"parent_id":40,"node_name":"LABORS-LION-H2-C1","updated_at":"2024-02-19 22:39:17 +0000","file_hash":None,"archive_id":None},
        ])

        # we have a config
        config = {
            "project_name": "labors",
            "config": {
                "bucket_name": "pics",
                "models": {
                    "victim": {
                        "scripts": [
                            {
                                "type": "file_collection",
                                "folder_path": "family",
                                "file_match": "{LABORS-LION-H2-C1,LABORS-LION-H2-C2}/family_photo.*.png",
                                "attribute_name": "family_photos"
                            }
                        ]
                    }
                }
            }
        }
        rules={
            'victim': r'^LABORS-LION-H\d+-C\d+$'
        }

        # given these we get an update
        update = MetisLoaderConfig(**config, rules=rules).update_for(tail)

        assert asdict(update)['revisions'] == {
            'victim': {
                'LABORS-LION-H2-C1': {
                    'family_photos': [
                        {
                            'original_filename': 'family_photo.1.png',
                            'path': 'metis://labors/pics/family/LABORS-LION-H2-C1/family_photo.1.png'
                        },
                        {
                            'original_filename': 'family_photo.2.png',
                            'path': 'metis://labors/pics/family/LABORS-LION-H2-C1/family_photo.2.png'
                        }
                    ]
                }
            }
        }

    def test_matches_tails_to_data_frame_scripts(self):
        # we have a tail
        tail = prep_tail('labors', 'pics', [
            {"type":"file","id":51,"parent_id":41,"node_name":"village-1.tsv","updated_at":"2023-11-11 22:39:17 +0000","file_hash":"8277e0910d750195b448797616e091ad","archive_id":None},
            {"type":"file","id":52,"parent_id":41,"node_name":"village-2.csv","updated_at":"2023-11-11 22:39:17 +0000","file_hash":"9277e0910d750195b448797616e091ad","archive_id":None},
            {"type":"parent","id":41,"parent_id":None,"node_name":"villages","updated_at":"2023-08-03 22:39:17 +0000","file_hash":None,"archive_id":None}
        ])

        # we have configs:
        # tsv version, file columns do match the map
        config1 = {
            "project_name": "labors",
            "config": {
                "bucket_name": "pics",
                "models": {
                    "victim": {
                        "scripts": [
                            {
                                "type": "data_frame",
                                "folder_path": "villages",
                                "file_match": "*.tsv",
                                "format": "tsv",
                                "column_map": {
                                    "name": "name",
                                    "species": "SPECIES"
                                }
                            }
                        ]
                    }
                }
            }
        }
        # csv format, file columns do match the map, script also set to test auto detection that its a csv file.
        config2 = {
            "project_name": "labors",
            "config": {
                "bucket_name": "pics",
                "models": {
                    "victim": {
                        "scripts": [
                            {
                                "type": "data_frame",
                                "folder_path": "villages",
                                "file_match": "*.csv",
                                "format": "auto-detect",
                                "column_map": {
                                    "name": "name",
                                    "species": "SPECIES"
                                }
                            }
                        ]
                    }
                }
            }
        }
        # tsv format, but column_map will not matching the file columns (Error check)
        config3 = {
            "project_name": "labors",
            "config": {
                "bucket_name": "pics",
                "models": {
                    "victim": {
                        "scripts": [
                            {
                                "type": "data_frame",
                                "folder_path": "villages",
                                "file_match": "*.tsv",
                                "format": "auto-detect",
                                "column_map": {
                                    "name": "name",
                                    "species": "species"
                                }
                            }
                        ]
                    }
                }
            }
        }
        # csv format, but tsv file, otherwise valid
        config4 = {
            "project_name": "labors",
            "config": {
                "bucket_name": "pics",
                "models": {
                    "victim": {
                        "scripts": [
                            {
                                "type": "data_frame",
                                "folder_path": "villages",
                                "file_match": "*.tsv",
                                "format": "csv",
                                "column_map": {
                                    "name": "name",
                                    "species": "SPECIES"
                                }
                            }
                        ]
                    }
                }
            }
        }
        # we have this rule for the target model' identifiers
        rules={
            'victim': r'^LABORS-LION-H\d+-C\d+$'
        }
        # we have this template for the target model
        models={
            'victim': Template(
                identifier='name',
                attributes={
                    'name': Attribute({ 'attribute_type': 'identifier' }),
                    'species': Attribute({ 'attribute_type': 'string' })
                }
            )
        }

        metis = mock.Mock()

        # And we have these files
        files = {
            'village-1.tsv': 'name\tSPECIES\nLABORS-LION-H2-C1\tlion\n',
            'village-2.csv': 'name,SPECIES\nLABORS-LION-H2-C1,lion\n'
        }

        def dummy_file(file):
            return StringIO(files[file.file_name])

        metis.open_file = mock.Mock(side_effect=dummy_file)

        # given these we get an update
        update1 = MetisLoaderConfig(**config1, rules=rules).update_for(tail, metis, models)
        update2 = MetisLoaderConfig(**config2, rules=rules).update_for(tail, metis, models)

        assert asdict(update1)['revisions'] == {
            'victim': {
                'LABORS-LION-H2-C1': {
                    'species': 'lion'
                }
            }
        }
        assert asdict(update1)['revisions'] == asdict(update2)['revisions']

        # Note that path3 also makes use of auto file format detection, but is expected to fail later, after the data_frame is read in and found to be missing a mapped column.
        with raises(MetisLoaderError, match=r"missing column.*species"):
            MetisLoaderConfig(**config3, rules=rules).update_for(tail, metis, models)

        with raises(MetisLoaderError, match=r"Check the 'format' configuration"):
            MetisLoaderConfig(**config4, rules=rules).update_for(tail, metis, models)

        # Now adjusting the model template to where the identifier would be missing from the column_map
        models_id={
            'victim': Template(
                identifier='id',
                attributes={
                    'name': Attribute({ 'attribute_type': 'identifier' }),
                    'species': Attribute({ 'attribute_type': 'string' })
                }
            )
        }
        with raises(MetisLoaderError, match=r"Identifier attribute is missing"):
            MetisLoaderConfig(**config1, rules=rules).update_for(tail, metis, models_id)

        # Now adjusting the model template to where column_map targets non-existent attributes
        models_missing={
            'victim': Template(
                identifier='name',
                attributes={
                    'name': Attribute({ 'attribute_type': 'identifier' })
                }
            )
        }
        with raises(MetisLoaderError, match=r"attribute\(s\) that don't exist: species"):
            MetisLoaderConfig(**config1, rules=rules).update_for(tail, metis, models_missing)

    def test_universal_linker_dag(self):
        pass
