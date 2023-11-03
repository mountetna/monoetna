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
        # For table version, DO blank past values
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
                                "blank_table": True,
                                "column_map": {
                                    "name": "name",
                                    "species": "SPECIES",
                                    "target_name": "target_name"
                                }
                            }
                        ]
                    }
                }
            }
        }
        # Same as above, just additionally set to ignore data that are "__" or empty strings
        config1hole = {
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
                                "blank_table": True,
                                "column_map": {
                                    "name": "name",
                                    "species": "SPECIES",
                                    "target_name": "target_name"
                                },
                                "values_to_ignore": "__,"
                            }
                        ]
                    }
                }
            }
        }
        # csv format, file columns do match the map, script also set to test auto detection that its a csv file.
        # For table version, DON'T blank past values
        # Config lacks a 'values_to_ignore' definition.
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
                                "blank_table": False,
                                "column_map": {
                                    "name": "name",
                                    "species": "SPECIES",
                                    "target_name": "target_name"
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
                                "blank_table": True,
                                "column_map": {
                                    "name": "name",
                                    "species": "species",
                                    "target_name": "target_name"
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
                                "blank_table": True,
                                "column_map": {
                                    "name": "name",
                                    "species": "SPECIES",
                                    "target_name": "target_name"
                                }
                            }
                        ]
                    }
                }
            }
        }
        # we have this rule for the target model' identifiers
        rules={
            'name': r'^LABORS-LION-H\d+-C\d+$',
            'victim': r'^LABORS-LION-H\d+-C\d+$'
        }
        # we have these potential template summaries for the target model
        model_std={
            'victim': {
                'template': Template(
                    identifier='name',
                    parent='stub',
                    attributes={
                        'name': Attribute({ 'attribute_type': 'identifier' }),
                        'species': Attribute({ 'attribute_type': 'string' }),
                        'target_name': Attribute({ 'attribute_type': 'string' })
                    }
                ),
                'isTable': False
            }
        }
        model_table={
            'victim': {
                'template': Template(
                    identifier='id',
                    parent='name',
                    attributes={
                        'name': Attribute({ 'attribute_type': 'parent' }),
                        'species': Attribute({ 'attribute_type': 'string' }),
                        'target_name': Attribute({ 'attribute_type': 'string' })
                    }
                ),
                'isTable': True
            }
        }

        metis = mock.Mock()

        # And we have these files
        files = {
            'village-1.tsv': 'name\tSPECIES\ttarget_name\nLABORS-LION-H2-C1\t__\t\n',
            'village-2.csv': 'name,SPECIES,target_name\nLABORS-LION-H2-C1,__,\n'
        }

        def dummy_file(file):
            return StringIO(files[file.file_name])

        metis.open_file = mock.Mock(side_effect=dummy_file)

        # given the 'good' configs, we get an update
        update1 = MetisLoaderConfig(**config1, rules=rules).update_for(tail, metis, model_std)
        update2 = MetisLoaderConfig(**config2, rules=rules).update_for(tail, metis, model_std)
        update3 = MetisLoaderConfig(**config1, rules=rules).update_for(tail, metis, model_table)
        update4 = MetisLoaderConfig(**config2, rules=rules).update_for(tail, metis, model_table)

        assert asdict(update1)['revisions'] == {
            'victim': {
                'LABORS-LION-H2-C1': {
                    'species': '__',
                    'target_name': ''
                }
            }
        }
        assert asdict(update1)['revisions'] == asdict(update2)['revisions']
        assert asdict(update3)['revisions'] == {
            'name': {
                'LABORS-LION-H2-C1': {'victim': ['::temp-id-0']}
            },
            'victim': {
                '::temp-id-0': {
                    'name': 'LABORS-LION-H2-C1',
                    'species': '__',
                    'target_name': ''
                }
            }
        }
        assert asdict(update4)['revisions'] == {
            'victim': {
                '::temp-id-0': {
                    'name': 'LABORS-LION-H2-C1',
                    'species': '__',
                    'target_name': ''
                }
            }
        }

        # given a 'good' config, with values_to_ignore matching a value in the file, we get a shortened update
        update5 = MetisLoaderConfig(**config1hole, rules=rules).update_for(tail, metis, model_std)
        assert asdict(update5)['revisions'] == {
            'victim': {
                'LABORS-LION-H2-C1': {}
            }
        }

        # Note that path3 also makes use of auto file format detection, but is expected to fail later, after the data_frame is read in and found to be missing a mapped column.
        with raises(MetisLoaderError, match=r"missing column.*species"):
            MetisLoaderConfig(**config3, rules=rules).update_for(tail, metis, model_std)
        with raises(MetisLoaderError, match=r"missing column.*species"):
            MetisLoaderConfig(**config3, rules=rules).update_for(tail, metis, model_table)

        with raises(MetisLoaderError, match=r"Check the 'format' configuration"):
            MetisLoaderConfig(**config4, rules=rules).update_for(tail, metis, model_std)
        with raises(MetisLoaderError, match=r"Check the 'format' configuration"):
            MetisLoaderConfig(**config4, rules=rules).update_for(tail, metis, model_table)

        # Now adjusting the model template to where the identifier would be missing from the column_map
        def models_id(isTable):
            return {
                'victim': {
                    'template': Template(
                        identifier='id',
                        attributes={
                            'name': Attribute({ 'attribute_type': 'identifier' }),
                            'species': Attribute({ 'attribute_type': 'string' }),
                            'target_name': Attribute({ 'attribute_type': 'string' })
                        }
                    ),
                    'isTable': isTable
                }
            }
        with raises(MetisLoaderError, match=r"Identifier attribute is missing"):
            MetisLoaderConfig(**config1, rules=rules).update_for(tail, metis, models_id(False))
        with raises(MetisLoaderError, match=r"Parent attribute is missing"):
            MetisLoaderConfig(**config1, rules=rules).update_for(tail, metis, models_id(True))

        # Now adjusting the model template to where column_map targets non-existent attributes
        def models_missing(isTable):
            return {
                'victim': {
                    'template': Template(
                        identifier='name',
                        attributes={
                            'name': Attribute({ 'attribute_type': 'identifier' }),
                            'target_name': Attribute({ 'attribute_type': 'string' })
                        }
                    ),
                    'isTable': False
                }
            }
        with raises(MetisLoaderError, match=r"attribute\(s\) that don't exist: species"):
            MetisLoaderConfig(**config1, rules=rules).update_for(tail, metis, models_missing(False))
        with raises(MetisLoaderError, match=r"attribute\(s\) that don't exist: species"):
            MetisLoaderConfig(**config1, rules=rules).update_for(tail, metis, models_missing(True))

    def test_universal_linker_dag(self):
        pass
