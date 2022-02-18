from airflow.utils.task_group import TaskGroup

from etna import (
    metis_etl,
    link,
    MetisEtlHelpers,
)
import re

timepoint_regex = re.compile(r'^(?P<timepoint>MVIR1-HS\d+-DN?\d+)[A-Z]+.*')
patient_regex = re.compile(r'^(?P<patient>MVIR1-HS\d+)-D.*')
day_regex = re.compile(r'.*-DN?(?P<day>\d+)$')
raw_bulk_rna_seq_record_regex = re.compile(r'^bulk_RNASeq/raw/[^/]*')
processed_bulk_rna_seq_record_regex = re.compile(r'^bulk_RNASeq/processed/[^/]*')

def link_rna_seq_with_timepoint(matches):
    for match in matches:
        timepoint_id = timepoint_regex.match(match.record_name)
        if not timepoint_id:
            continue
        timepoint_id = timepoint_id.group('timepoint')

        patient_id = patient_regex.match(match.record_name)
        if not patient_id:
            continue
        patient_id = patient_id.group('patient')

        day = day_regex.match(timepoint_id)
        if not day:
            continue

        day = int(day.group('day'))

        if 'DN' in timepoint_id:
            day *= -1

        update = match.as_update(dict(timepoint=timepoint_id))
        update.update_record('timepoint', timepoint_id, dict(patient=patient_id, day=day))
        yield match, update


@metis_etl('mvir1', 'data', version=13)
def mvir1_data_metis_etl(helpers: MetisEtlHelpers):
    with TaskGroup("bulk_rna_seq_raw"):
        matches = helpers.find_record_folders('rna_seq', raw_bulk_rna_seq_record_regex)
        matches = link(dry_run=True)(link_rna_seq_with_timepoint)(matches)
        listed_matches = helpers.list_match_folders(matches)

        helpers.link_matching_files(listed_matches, 'raw_fastq_files', file_regex=re.compile(r'.*\.fastq\.gz$'), dry_run=True)

    with TaskGroup("bulk_rna_seq_processed"):
        matches = helpers.find_record_folders('rna_seq', processed_bulk_rna_seq_record_regex)
        matches = helpers.filter_by_timur(matches)
        listed_matches = helpers.list_match_folders(matches)

        # helpers.link_matching_file(listed_matches, )
