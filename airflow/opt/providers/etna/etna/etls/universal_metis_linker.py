from etna.hooks.etna import EtnaHook

class UniversalMetisLinker:
    """
    Creates a universal linker etl. This consumes tails on metis, but it does
    it for a set of projects interactively
    """

    def __init__(self):
        self.hook = EtnaHook.for_project('administration')
        self.hook.use_token_auth = True

    def run(self):
        # Linker steps:
        # 1. Get ETL configs from Polyphemus for loader 'metis'.
        hook = self.hook

        with hook.polyphemus() as polyph:
            etls = polyph.etl_list(job_type='metis')

        project_names = set(etl.project_name for etl in running_etls)

        # 2. Get Metis tail for the specified time period and the projects to run from ETL configs

        project_revisions = {}

        for etl in etls:
            if not should_run(etl):
                continue

            project_revisions.setdefault( etl.project_name, {} )
            try:
                for model_name, model in etl.config.items():
                    for script in model['scripts']:
                        tail = ensure_tail(tail_collection)

                        if script['type'] == 'file':
                            files = find_tail_files(
                                tail,
                                script['folder_path'],
                                script['file_match']
                            )
                            for file in files:
                                record_name = get_identifier( model_name, file )
                                add_to_revision(
                                    project_revisions[ etl.project_name ],
                                    model_name,
                                    record_name,
                                    script['attribute_name'],
                                    file_reference_for(file)
                                )
                        elif script['type'] == 'file_collection':
                            files = find_tail_files(
                                tail,
                                script['folder_path'],
                                script['file_match']
                            )

                            # group all files by parent folders
                            folders = group_by_parent_folder(files)
                            for folder in folders:
                                record_name = get_identifier( model_name, file )
                                add_to_revision(
                                    project_revisions[ etl.project_name ],
                                    model_name,
                                    record_name,
                                    script['attribute_name'],
                                    [
                                        file_reference_for(file)
                                        for file in folder
                                    ]
                                )
                        elif script['type'] == 'data_frame':
                            files = find_tail_files()
            except:
                pass

        with hook.metis() as metis:
            for project_name, bucket_name in project_buckets:
                return pickled(load_metis_files_batch(metis, bucket_name))

        # 3. Match items from ETL configs in the tail

        # 4. Compose updating using matched items and send to Magma
