from airflow.decorators import dag, task
from datetime import datetime, timedelta

from etna.hooks.etna import EtnaHook
from airflow.operators.python import get_current_context
from etna.etls.etl_task_batching import get_batch_range

from mountetna import UpdateRequest
from mountetna.polyphemus import EtlConfigResponse
from dataclasses import dataclass, field, asdict
from typing import Dict, Any
from serde import serialize, deserialize
import re
import itertools
import pandas
import numpy
import secrets

import logging

logger = logging.getLogger(__name__)

class MetisLoaderError(Exception):
    pass

STAR_MAGIC = secrets.token_hex(12)
BRACE1_MAGIC = secrets.token_hex(12)
BRACE2_MAGIC = secrets.token_hex(12)

@serialize
@deserialize
@dataclass
class MetisLoaderConfig(EtlConfigResponse):
    rules: Dict[ str, str ] = field(default_factory=dict)

    @property
    def bucket_key(self):
        return (self.project_name, self.config['bucket_name'])

    def script_files(self, script, tail):
        try:
            # simulate glob
            match = script['folder_path'] + '/**/' + script['file_match']
            match = re.sub(r'\*', STAR_MAGIC, match)
            match = re.sub(r'\{', BRACE1_MAGIC, match)
            match = re.sub(r'\}', BRACE2_MAGIC, match)
            match = re.escape(match)

            # **/ matches anything including nothing
            match = re.sub(STAR_MAGIC*2+'/', '(?:([^/]*(/|$))*)', match)
      
            # /** matches any number of directories
            match = re.sub(f'/{STAR_MAGIC*2}$', '(?:(/[^/]*)*)', match)

            # * matches any non-slash including nothing
            match = re.sub(STAR_MAGIC, '[^/]*', match)

            # {a,b,c,d} matches a group
            match = re.sub(f"{BRACE1_MAGIC}((?:[^,]+,)+[^,]+){BRACE2_MAGIC}", lambda m : f"({m[1].replace(',','|')})", match)
            match = re.compile(match)
        except:
            raise MetisLoaderError(f"file_match is not a valid regular expression: {script['file_match']}.")

        files = [
            file
            for file in tail[0]
            if match.match(file.file_path)
        ]
        return files

    def data_frame_update(self, model_name, script, tail, update, metis, model):
        isTable=model['isTable']
        template=model['template']
        model_attributes = list(template.attributes.keys())
        id = template.identifier

        column_map = script['column_map']
        attributes = list(column_map.keys())
        columns = list(column_map.values())

        if not isTable and not id in attributes:
            raise MetisLoaderError(f"Identifier attribute is missing from the 'column_map' of {model_name} data_frame loader.")
        if isTable and not template.parent in attributes:
            raise MetisLoaderError(f"Parent attribute is missing from the 'column_map' of {model_name} data_frame loader.")
        if not set(attributes).issubset(model_attributes):
            missing = ', '.join([attr for attr in attributes if attr not in model_attributes])
            raise MetisLoaderError(f"'column_map' of {model_name} data_frame loader targets attribute(s) that don't exist: {missing}.")

        # Parse Files
        files = self.script_files(script, tail)
        for file in files:
            data = None
            with metis.open_file(file) as file_reader:
                if script['format']=="csv":
                    separator = ","
                if script['format']=="tsv":
                    separator = "\t"
                if script['format']=="auto-detect":
                    if re.search("csv$", file.file_path) is not None:
                        separator = ","
                    elif re.search("tsv$", file.file_path) is not None:
                        separator = "\t"
                    else:
                        # Attempt to use pandas' automated detection
                        separator = None
                # Load dataframe
                data = pandas.read_table(file_reader, sep=separator, engine = 'python', na_filter=False)
            if len(data.columns) < 2:
                raise MetisLoaderError(f"{file.file_name} seems to have fewer than 2 columns. Check the 'format' configuration for this data_frame loader.")
            if not set(columns).issubset(data.columns):
                missing = ', '.join([col for col in columns if col not in data.columns])
                raise MetisLoaderError(f"{file.file_name} is missing column(s) targetted by {model_name} data_frame loader 'column_map': {missing}.")
            if pandas.isna(data).values.any():
                raise MetisLoaderError(f"{file.file_name} has unexpected NA values after all parsing. Data rows may be shorter than the column row indicates.")
            # Trim to mapped columns and convert to attribute names
            data = data.rename(columns={v: k for k,v in column_map.items()})[attributes]
            # Blank data equaling values_to_ignore by setting as None (None's are ignored when updates are constructed)
            if 'values_to_ignore' in script:
                replacements = {k: None for k in script['values_to_ignore'].split(',')}
                data=data.replace(replacements)
            # Determine Updates
            if isTable:
                data['__temp__']=['::temp-id-' + str(temp) for temp in data.index]
                data = data.set_index('__temp__', drop=True)
                if script.get('blank_table', False):
                    for parent_name in set(list(data[template.parent])):
                        if not self.get_identifier(template.parent, parent_name):
                            continue
                        update.update_record(template.parent, parent_name, {model_name: list(data[data[template.parent]==parent_name].index)})
                for name, attributes in data.T.to_dict().items():
                    update.update_record(model_name, name, {k: v for k,v in attributes.items() if v is not None})
            else:
                data = data.set_index(id, drop=True)
                for name, attributes in data.T.to_dict().items():
                    if model_name in self.rules and not self.get_identifier(model_name, str(name)):
                        continue
                    update.update_record(model_name, name, {k: v for k,v in attributes.items() if v is not None})

    def file_collection_update(self, model_name, script, tail, update):
        if model_name not in self.rules:
            raise f"Cannot filter by file without a rule for {model_name}"

        files = self.script_files(script, tail)

        named_files = itertools.groupby(
            files,
            lambda f : self.get_identifier(model_name, f.file_path)
        )
        for name, name_files in named_files:
            if not name:
                continue
            update.update_record(
                model_name,
                name,
                {
                    script['attribute_name']: [
                        file.as_magma_file_attribute
                        for file in name_files
                    ]
                }
            )


    def get_identifier(self, model_name, file_name):
        match = re.search(self.rules[model_name][1:-1], file_name)

        if not match: 
            return None

        return match.group(0)

    def file_update(self, model_name, script, tail, update):
        if model_name not in self.rules:
            raise f"Cannot filter by file without a rule for {model_name}"

        files = self.script_files(script, tail)
        for file in files:
            name = self.get_identifier( model_name, file.file_path )
            if not name:
                continue
            update.update_record(
                model_name,
                name,
                {
                    script['attribute_name']: file.as_magma_file_attribute
                }
            )

    def update_for(self, tail, metis=None, models=None):
        update = UpdateRequest(project_name=self.project_name, dry_run=(not self.params.get('commit', False)), autolink = self.config.get('autolink', False))

        for model_name, model_config in self.config.get("models",{}).items():
            for script in model_config.get("scripts",[]):
                if script['type'] == 'file':
                    self.file_update( model_name, script, tail, update )
                elif script['type'] == 'file_collection':
                    self.file_collection_update( model_name, script, tail, update )
                elif script['type'] == 'data_frame':
                    self.data_frame_update( model_name, script, tail, update, metis, model=models[model_name])
                else:
                    raise MetisLoaderError(f"Invalid type for script {script['type']} for model {model_name}")

        return update

def MetisLinker(interval_minutes=5):
    @dag(
        max_active_runs=1,
        schedule_interval=timedelta(minutes=interval_minutes),
        start_date=datetime(2023, 4, 12),
        default_args={
            'depends_on_past': False,
            'wait_for_downstream': False
        }
    )
    def MetisLinker():
        hook = EtnaHook('etna_administration', use_token_auth=True)
        
        @task
        def read_configs():
            # read configs from polyphemus
            context = get_current_context()
            start, end = get_batch_range(context)

            with hook.polyphemus() as polyphemus:
                configs_list = polyphemus.list_all_etl_configs(job_type='metis')
                configs_list = [ MetisLoaderConfig(**asdict(config)) for config in configs_list.configs if config.should_run(start,end) ]

                for config in configs_list:
                    polyphemus.etl_update(
                        project_name=config.project_name,
                        config_id=config.config_id,
                        status='running',
                        ran_at=datetime.now().isoformat()
                    )
                return configs_list
        
        @task
        def get_rules(configs):
            with hook.gnomon() as gnomon:
                project_names = [ config.project_name for config in configs ]
                rules_list = gnomon.rules(project_names=project_names)
                return rules_list.rules

        @task
        def get_models(configs):
            def get_template(project, model, magma, check_table=True):
                raw_return=magma.retrieve(project_name=project, model_name=model, record_names=[], attribute_names="all", hide_templates=False)
                raw_template=raw_return.models[model].template
                output = {
                    'template': raw_template,
                    'isTable': False
                    }
                # Add isTable
                if check_table:
                    # Shove as 'IsTable' attribute (not snake case so should never conflict!)
                    output['isTable'] = raw_template.parent!='project' and get_template(project, raw_template.parent, magma, False)['template'].attributes[model].attribute_type=='table'
                return output
            # Minimize template grabs because some templates are LARGE
            project_names = set([ config.project_name for config in configs ])
            project_models = {p: [] for p in project_names}
            for config in configs:
                project_models[config.project_name].extend(list(config.config.get("models",{}).keys()))
            with hook.magma() as magma:
                return {
                    project: {model: get_template(project, model, magma) for model in list(set(project_models[project]))}
                    for project in project_names
                }

        @task
        def collect_tails(configs):
            # get project/bucket pairs from configs
            context = get_current_context()
            start, end = get_batch_range(context)
            print(f"Collecting tails in range {str(dict(start=start, end=end))}")

            buckets = dict()
            for config in configs:
                ran_at = datetime.fromisoformat(config.ran_at or '2010-01-01')
                if config.bucket_key not in buckets or buckets[ config.bucket_key ] > ran_at:
                    buckets[ config.bucket_key ] = ran_at
                    
            with hook.metis() as metis:
                tails = dict()
                for (project_name, bucket_name), ran_at in buckets.items():
                    try:
                        files = metis.tail(
                            project_name=project_name,
                            bucket_name=bucket_name,
                            type='files',
                            batch_start=ran_at,
                            batch_end=end
                        )
                        tails[ (project_name, bucket_name) ] = { 'files': files }
                    except Exception as error:
                        logging.error(error, exc_info=True)
                        tails[ (project_name, bucket_name) ] = { 'error': repr(error) }
                return tails

        @task
        def process_tails(configs, tails, rules, models):
            updates = {}
            for config in configs:
                try:
                    # find associated tails
                    tail = tails[ config.bucket_key ]
                    
                    if 'error' in tail:
                        raise Exception(tail['error'])

                    if config.project_name not in rules:
                        raise MetisLoaderError(f"No rules found for {config.project_name}")

                    config.rules = rules[ config.project_name ]

                    with hook.metis() as metis:
                        updates[ config.config_id ] = config.update_for(tail['files'], metis, models[config.project_name])
                except Exception as error:
                    logging.error(error, exc_info=True)
                    updates[ config.config_id ] = { 'error': repr(error) }

            return updates

        def error_report(**i):
            summary = f'''
===============================
Error : {i['start']} -> {i['end']} 
Job failed with error: {i['error']}
==============================='''

            return summary

        def upload_report(**i):
            summary = f'''
===============================
Upload Summary : {i['start']} -> {i['end']} 
Models: {', '.join(i['response'].models.keys())}
Committed to Magma: {i['commit']}
Autolinked Parent Identifiers: {i['autolink']}
'''
            for model_name, model in i['response'].models.items():
                summary += f"{model_name} records updated: {', '.join(model.documents.keys())}\n"

            summary += "==============================="

            return summary
                    
        @task
        def post_updates(configs, updates):
            context = get_current_context()
            start, end = get_batch_range(context)
            for config in configs:
                try:
                    if 'error' in updates[ config.config_id ]:
                        raise Exception(updates[ config.config_id ]['error'])

                    with hook.magma() as magma:
                        response = magma.update(updates[ config.config_id ], page_size=1000)

                    with hook.polyphemus() as polyphemus:
                        state = {}
                        if config.run_interval == 0:
                            state['run_interval'] = -1
                            state['params'] = {}

                        polyphemus.etl_update(
                            project_name=config.project_name,
                            config_id=config.config_id,
                            status='completed',
                            **state
                        )
                        polyphemus.etl_add_output(
                            project_name=config.project_name,
                            config_id=config.config_id,
                            append=True,
                            output=upload_report(
                                project_name=config.project_name,
                                response=response,
                                commit=config.params.get('commit', False),
                                autolink=config.config.get('autolink', False),
                                start=start,
                                end=end
                            )
                        )
                except Exception as error:
                    logging.error(error, exc_info=True)
                    with hook.polyphemus() as polyphemus:
                        polyphemus.etl_update(
                            project_name=config.project_name,
                            config_id=config.config_id,
                            status='error',
                            run_interval=-1
                        )
                        polyphemus.etl_add_output(
                            project_name=config.project_name,
                            config_id=config.config_id,
                            append=True,
                            output=error_report(start=start,end=end,error=repr(error))
                        )
            
        configs = read_configs()
        tails = collect_tails(configs)
        rules = get_rules(configs)
        models = get_models(configs)
        updates = process_tails(configs, tails, rules, models)
        post_updates(configs, updates)
    return MetisLinker()
