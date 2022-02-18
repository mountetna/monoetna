# from datetime import timedelta
# from etna import (
#     run_on_docker,
#     metis_etl,
#     filter_by_record_directory,
#     filter_by_exists_in_timur,
#     list_contents_of_matches,
#     pickled,
#     link,
# )
# from airflow.decorators import task
# from datetime import datetime
# import re
#
# timepoint_regex = re.compile(r'^(?P<timepoint>MVIR1-HS\d+-DN?\d+)[A-Z]+.*')
# patient_regex = re.compile(r'^(?P<patient>MVIR1-HS\d+)-D.*')
# day_regex = re.compile(r'.*-DN?(?P<day>\d+)$')
# bulk_rna_seq_record_regex = re.compile(r'^bulk_RNASeq/raw/[^/]*')
#
# pool_identifier_tri_regex = re.compile(r'^([^.]*\.[^.]*\.)(.*)$')
# non_pool_identifier_tri_regex = re.compile('^([^_.]*[_.][^_.]*[_.])([^_.]*)(.*)$')
#
#
# def correct_ipi_record_name(name):
#     if 'POOL' in name:
#         name = name.replace('_', '.')
#         match = pool_identifier_tri_regex.match(name)
#         if match:
#             name = match.groups(1) + match.groups(2).lower()
#     else:
#         match = non_pool_identifier_tri_regex.match(name)
#         if match:
#             name = match.groups(1).replace('_', '.') + match.groups(2).lower() + match.groups(3)
#     return name
#
#
# assert correct_ipi_record_name('IPIPOOL001_P1_scRNA_XVIPBMC') == "IPIPOOL001.P1.scrna.xvipbmc"
# assert correct_ipi_record_name('IPILUNG094.N1.scrna.CD45_enriched') == "IPILUNG094.N1.scrna.CD45_enriched"
# assert correct_ipi_record_name('IPIPOOL001_P1_scRNA_XVIPBMC') == 'IPIPOOL001.P1.scrna.xvipbmc'
#
# single_cell_root = r'^single_cell_[^/]*'
# # Pools are nested
# pool_root = r'.*/IPIPOOL[^/]+/IPIPOOL[^/]+'
# non_pool_root = r'.*/IPI((?!POOL)[^/])+'
#
#
# def _fastq_matches(listed_matches):
#     for match, files in listed_matches:
#         fastqs = [f for f in files if f.file_path.endswith(".fastq.gz")]
#         # Do not 'undo' fastq linkage if we happen to, say, match an empty folder.
#         if fastqs:
#             yield match, fastqs
#
# # Initializes tasks needed for processing the pooled single cell fastq files
# def pooled_sc_fastq_tasks(tail_folders, tail_files, hook):
#     @task
#     def find_pooled_sc_raw_folders(folders, files):
#         matches = filter_by_record_directory(folders + files, re.compile(f"{single_cell_root}/raw/{pool_root}"),
#                                              corrected_record_name=correct_ipi_record_name)
#         with hook.magma() as magma:
#             matches = filter_by_exists_in_timur(magma, matches, 'sc_rna_seq_pool')
#         with hook.metis() as metis:
#             return pickled(list_contents_of_matches(metis, matches))
#
#     @link('sc_rna_seq_pool', 'raw_fastq')
#     def link_sc_rna_seq_pooled_raw_fastqs(listed_matches):
#         return _fastq_matches(listed_matches)
#
#     link_sc_rna_seq_pooled_raw_fastqs(find_pooled_sc_raw_folders(tail_folders, tail_files))
#
# # Initializes tasks needed for processing the non pooled single cell fastq files
# def non_pooled_sc_fastq_tasks(tail_folders, tail_files, hook):
#     @task
#     def find_non_pooled_sc_raw_folders(folders, files):
#         matches = filter_by_record_directory(folders + files, re.compile(f"{single_cell_root}/raw/{non_pool_root}"),
#                                              corrected_record_name=correct_ipi_record_name)
#         with hook.magma() as magma:
#             matches = filter_by_exists_in_timur(magma, matches, 'sc_rna_seq')
#         with hook.metis() as metis:
#             return pickled(list_contents_of_matches(metis, matches))
#
#     @link('sc_rna_seq', 'raw_fastq')
#     def link_sc_rna_seq_non_pooled_raw_fastqs(listed_matches):
#         return _fastq_matches(listed_matches)
#
#     link_sc_rna_seq_non_pooled_raw_fastqs(find_non_pooled_sc_raw_folders(tail_folders, tail_files))
#
# def _sc_processed_folder_attributes(match):
#     pass
#
# def pooled_sc_processed_tasks(tail_folders, tail_files, hook):
#     @task
#     def find_pooled_sc_processed_folders(folders, files):
#         matches = filter_by_record_directory(folders + files, re.compile(f"{single_cell_root}/processed/{pool_root}"),
#                                              corrected_record_name=correct_ipi_record_name)
#         with hook.magma() as magma:
#             matches = filter_by_exists_in_timur(magma, matches, 'sc_rna_seq_pool')
#         with hook.metis() as metis:
#             return pickled(list_contents_of_matches(metis, matches))
#
#     @link('sc_rna_seq_pool', 'cite_antibody_key')
#     def link_sc_rna_seq_pooled_processed_files(listed_matches):
#         for match, files in listed_matches:
#             for file in files:
#                 if re.compile(r'^ADT_keep_features\.list$').match(file.file_name):
#                     yield match, file
#
#     @link('sc_rna_seq_pool', 'mux_index_key')
#     def link_sc_rna_seq_pooled_processed_files(listed_matches):
#         for match, files in listed_matches:
#             for file in files:
#                 if re.compile(r'^IDX_map\.tsx$').match(file.file_name):
#                     yield match, file
#
#     @link('sc_rna_seq_pool', 'processed_robject_rdata')
#     def link_sc_rna_seq_pooled_processed_files(listed_matches):
#         for match, files in listed_matches:
#             for file in files:
#                 if re.compile(r'.*_scTransformed_processed\.RData$').match(file.file_name):
#                     yield match, file
#
#     @link('sc_rna_seq_pool', 'processed_umap')
#     def link_sc_rna_seq_pooled_processed_files(listed_matches):
#         for match, files in listed_matches:
#             for file in files:
#                 if re.compile(r'.*_umap\.pdf$').match(file.file_name):
#                     yield match, file
#
#     @link('sc_rna_seq_pool', 'processing_pipeline_parameters')
#     def link_sc_rna_seq_pooled_processed_files(listed_matches):
#         for match, files in listed_matches:
#             for file in files:
#                 if re.compile(r'^cutoffs\.yml$').match(file.file_name):
#                     yield match, file
#
#     @link('sc_rna_seq_pool', 'filtered_counts_h5')
#     def link_sc_rna_seq_pooled_processed_files(listed_matches):
#         for match, files in listed_matches:
#             for file in files:
#                 if re.compile(r'^filtered_feature_bc_matrix\.h5$').match(file.file_name):
#                     yield match, file
#
#     @link('sc_rna_seq_pool', 'raw_counts_h5')
#     def link_sc_rna_seq_pooled_processed_files(listed_matches):
#         for match, files in listed_matches:
#             for file in files:
#                 if re.compile(r'^raw_feature_bc_matrix\.h5$').match(file.file_name):
#                     yield match, file
#
#     @link('sc_rna_seq_pool', 'tenx_aligned_bam')
#     def link_sc_rna_seq_pooled_processed_files(listed_matches):
#         for match, files in listed_matches:
#             for file in files:
#                 if re.compile(r'^possorted_genome_bam\.bam$').match(file.file_name):
#                     yield match, file
#
#     @link('sc_rna_seq_pool', 'tenx_aligned_bam_index')
#     def link_sc_rna_seq_pooled_processed_files(listed_matches):
#         for match, files in listed_matches:
#             for file in files:
#                 if re.compile(r'^possorted_genome_bam\.bam\.bai$').match(file.file_name):
#                     yield match, file
#
#     @link('sc_rna_seq_pool', 'tenx_cloupe_file')
#     def link_sc_rna_seq_pooled_processed_files(listed_matches):
#         for match, files in listed_matches:
#             for file in files:
#                 if re.compile(r'^cloupe\.cloupe$').match(file.file_name):
#                     yield match, file
#
#     @link('sc_rna_seq_pool', 'tenx_metrics_csv')
#     def link_sc_rna_seq_pooled_processed_files(listed_matches):
#         for match, files in listed_matches:
#             for file in files:
#                 if re.compile(r'^metrics_summary\.csv$').match(file.file_name):
#                     yield match, file
#
#     @link('sc_rna_seq_pool', 'tenx_molecule_info_h5')
#     def link_sc_rna_seq_pooled_processed_files(listed_matches):
#         for match, files in listed_matches:
#             for file in files:
#                 if re.compile(r'^molecule_info\.h5$').match(file.file_name):
#                     yield match, file
#
#     @link('sc_rna_seq_pool', 'tenx_web_summary')
#     def link_sc_rna_seq_pooled_processed_files(listed_matches):
#         for match, files in listed_matches:
#             for file in files:
#                 if re.compile(r'^web_summary\.html$').match(file.file_name):
#                     yield match, file
#
# @metis_etl('ipi', 'data', datetime(2022, 1, 1), 11)
# def ipi_data_metis_etl(tail_folders, tail_files, hook):
#     pooled_sc_fastq_tasks(tail_folders, tail_files, hook)
#     non_pooled_sc_fastq_tasks(tail_folders, tail_files, hook)