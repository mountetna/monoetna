import io

from airflow.utils.task_group import TaskGroup

from etna import (
    metis_etl,
    link,
    MetisEtlHelpers,
    UpdateRequest,
)
import re

from etna.etls.metis import MatchedRecordFolder
from etna.hooks.etna import File, EtnaHook, Metis

timepoint_regex = re.compile(r"^(?P<timepoint>MVIR1-HS\d+-DN?\d+)[A-Z]+.*")
patient_regex = re.compile(r"^(?P<patient>MVIR1-HS\d+)-D.*")
day_regex = re.compile(r".*-DN?(?P<day>\d+)$")
raw_bulk_rna_seq_record_regex = re.compile(r"^bulk_RNASeq/raw/[^/]*")
processed_bulk_rna_seq_record_regex = re.compile(r"^bulk_RNASeq/processed/[^/]*")


def link_rna_seq_with_timepoint(matches):
    for match in matches:
        timepoint_id = timepoint_regex.match(match.record_name)
        if not timepoint_id:
            continue
        timepoint_id = timepoint_id.group("timepoint")

        patient_id = patient_regex.match(match.record_name)
        if not patient_id:
            continue
        patient_id = patient_id.group("patient")

        day = day_regex.match(timepoint_id)
        if not day:
            continue

        day = int(day.group("day"))

        if "DN" in timepoint_id:
            day *= -1

        update = match.as_update(dict(timepoint=timepoint_id))
        update.update_record(
            "timepoint", timepoint_id, dict(patient=patient_id, day=day)
        )
        yield match, update


def process_gene_count_matrix(
    metis: Metis, match: MatchedRecordFolder, file: File
) -> UpdateRequest:
    with metis.open_file(file) as file:
        return match.as_update(produce_gene_matrix(file))


def produce_gene_matrix(gene_tab: io.BufferedReader):
    hook = EtnaHook.for_project()
    return {}


# produce_gene_matrix(io.BytesIO(b"asdfasdfadsf"))


def process_marked_dup_metrics(
    metis: Metis, match: MatchedRecordFolder, file: File
) -> UpdateRequest:
    with metis.open_file(file) as metis_file:
        return match.as_update(produce_marked_dup_metric(metis_file, file.file_path))


def produce_marked_dup_metric(marked_file: io.BufferedReader, file_path: str):
    line = marked_file.readline()
    while line:
        line = line.strip()
        if b"PERCENT_DUPLICATION" in line:
            dup_index = line.split(b"\t").index(b"PERCENT_DUPLICATION")
            break
        line = marked_file.readline()
    else:
        raise ValueError(
            f"Could not find PERCENT_DUPLICATION header in marked dup metrics file: {file_path}"
        )

    data_line = marked_file.readline().strip().split(b"\t")
    if len(data_line) <= dup_index:
        raise ValueError(
            f"Could not PERCENT_DUPLICATION value entry in marked dup metrics file: {file_path}"
        )

    return dict(duplication_pct=float(data_line[dup_index]))


# An example output tested
assert (
    produce_marked_dup_metric(
        io.BytesIO(
            b"""
## htsjdk.samtools.metrics.StringHeader
# MarkDuplicates INPUT=[Alignments/STAR211025_121627/415-1001-D4-RLT1out/Aligned.sortedByCoord.out.bam] OUTPUT=Alignments/STAR211025_121627/415-1001-D4-RLT1out/Aligned.sortedByCoord.duplicateRemoved.out.bam METRICS_FILE=Alignments/STAR211025_121627/415-1001-D4-RLT1out/marked_dup_metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true TMP_DIR=[Alignments/STAR211025_121627/415-1001-D4-RLT1out]    MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 SORTING_COLLECTION_SIZE_RATIO=0.25 TAG_DUPLICATE_SET_MEMBERS=false REMOVE_SEQUENCING_DUPLICATES=false TAGGING_POLICY=DontTag CLEAR_DT=true DUPLEX_UMI=false ADD_PG_TAG_TO_READS=true DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES PROGRAM_RECORD_ID=MarkDuplicates PROGRAM_GROUP_NAME=MarkDuplicates READ_NAME_REGEX=<optimized capture of last three ':' separated fields as numeric values> OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 MAX_OPTICAL_DUPLICATE_SET_SIZE=300000 VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false GA4GH_CLIENT_SECRETS=client_secrets.json USE_JDK_DEFLATER=false USE_JDK_INFLATER=false
## htsjdk.samtools.metrics.StringHeader
# Started on: Wed Nov 03 15:09:22 PDT 2021

## METRICS CLASS	picard.sam.DuplicationMetrics
LIBRARY	UNPAIRED_READS_EXAMINED	READ_PAIRS_EXAMINED	SECONDARY_OR_SUPPLEMENTARY_RDS	UNMAPPED_READS	UNPAIRED_READ_DUPLICATES	READ_PAIR_DUPLICATES	READ_PAIR_OPTICAL_DUPLICATES	PERCENT_DUPLICATION	ESTIMATED_LIBRARY_SIZE
Unknown Library	0	32863893	0	0	0	18525954	461915	0.563718	16764564

## HISTOGRAM	java.lang.Double
BIN	CoverageMult	all_sets	optical_sets	non_optical_sets
1.0	1.0046	11473929	0	11650466
2.0	1.146061	1919006	288694	1788982
"""
        ),
        "",
    )
    == dict(duplication_pct=float(b"0.563718"))
)


@metis_etl("mvir1", "data", version=13)
def mvir1_data_metis_etl(helpers: MetisEtlHelpers):
    with TaskGroup("bulk_rna_seq_raw"):
        matches = helpers.find_record_folders("rna_seq", raw_bulk_rna_seq_record_regex)
        matches = link(dry_run=True)(link_rna_seq_with_timepoint)(matches)
        listed_matches = helpers.list_match_folders(matches)

        helpers.link_matching_files(
            listed_matches,
            "raw_fastq_files",
            file_regex=re.compile(r".*\.fastq\.gz$"),
            dry_run=True,
        )

    with TaskGroup("bulk_rna_seq_processed"):
        matches = helpers.find_record_folders(
            "rna_seq", processed_bulk_rna_seq_record_regex
        )
        matches = helpers.filter_by_timur(matches)
        listed_matches = helpers.list_match_folders(matches)

        helpers.link_matching_file(
            listed_matches, "genome_alignment", re.compile(r".*\.bam$")
        )

        helpers.process_and_link_matching_file(
            listed_matches,
            re.compile(".*\.tab"),
            process_gene_count_matrix,
            dry_run=True,
        )
        helpers.process_and_link_matching_file(
            listed_matches,
            re.compile(".*marked_dup_metrics\.txt$"),
            process_marked_dup_metrics,
            dry_run=True,
        )
