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
            match = re.sub(STAR_MAGIC*2+'/', '(?:.*/|)', match)
            # * matches any non-slash including nothing
            match = re.sub(STAR_MAGIC, '[^/]*', match)
            # {a,b,c,d} matches a group
            match = re.sub(f"{BRACE1_MAGIC}((?:[^,]+,)+,[^,]+){BRACE2_MAGIC}", lambda m : f"({m[1].replace(',','|')})", match)
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
        files = self.script_files(script, tail)
        for file in files:
            with metis.open_file(file) as file_reader:
                data = pandas.read_csv(file_reader).replace({numpy.nan: None})
                for tube_name, table in data.groupby(model_name):
                    for row in table.to_dict('records'):
                        match = re.search(r'^(MVIR1-HS\d+)-[DM]N?[0-9]+', tube_name)
                        timepoint_name = match.group(0)
                        patient_name = match.group(1)
                        update.append_table('immunoassay', tube_name, 'analyte', row, 'analyte')
                        update.update_record('immunoassay', tube_name, { 'timepoint':  timepoint_name })
                        # These weird "M" records denote monthly timepoints and only exist in this CSV
                        # To be removed when we can do valid identifier interpolation
                        if re.search(r'-M[0-9]+$', timepoint_name):
                            update.update_record('timepoint', timepoint_name, { 'patient': patient_name })

    def file_collection_update(self, model_name, script, tail, update):
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
        match = re.search( self.rules[model_name][1:-1], file_name)

        if not match: 
            return None

        return match.group(0)

    def file_update(self, model_name, script, tail, update):
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

    def update_for(self, tail, metis=None):
        update = UpdateRequest(project_name=self.project_name, dry_run=(not self.params.get('commit', False)))

        for model_name, model_config in self.config.get("models",{}).items():
            for script in model_config.get("scripts",[]):
                if script['type'] == 'file':
                    self.file_update( model_name, script, tail, update )
                elif script['type'] == 'file_collection':
                    self.file_collection_update( model_name, script, tail, update )
                else:
                    raise MetisLoaderError(f"Invalid type for script {script['type']} for model {model_name}")

        return update

@dag(
    max_active_runs=1,
    schedule_interval=timedelta(minutes=60),
    start_date=datetime(2023, 4, 12),
    default_args={
        'depends_on_past': True,
        'wait_for_downstream': True
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
    def collect_tails(configs):
        # get project/bucket pairs from configs
        context = get_current_context()
        start, end = get_batch_range(context)
        print(f"Collecting tails in range {str(dict(start=start, end=end))}")

        buckets = [
            config.bucket_key
            for config in configs
            for model, model_config in config.config["models"].items()
            for script in model_config["scripts"]
        ]
                
        with hook.metis() as metis:
            tails = dict()
            for project_name, bucket_name in buckets:
                try:
                    files = metis.tail(
                        project_name=project_name,
                        bucket_name=bucket_name,
                        type='files',
                        batch_start=start,
                        batch_end=end
                    )
                    tails[ (project_name, bucket_name) ] = { 'files': files }
                except Exception as error:
                    tails[ (project_name, bucket_name) ] = { 'error': repr(error) }
            return tails

    @task
    def process_tails(configs, tails, rules):
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

                updates[ config.config_id ] = config.update_for(tail['files'])
            except Exception as error:
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
'''
        for model_name, model in i['response'].models.items():
            summary += f"{model_name} records updated: {len(model.documents)}\n"

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

                commit = config.params.get('commit', False)
                with hook.magma() as magma:
                    response = magma.update(updates[ config.config_id ])

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
                            commit=commit,
                            start=start,
                            end=end
                        )
                    )
            except Exception as error:
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
    updates = process_tails(configs, tails, rules)
    post_updates(configs, updates)
