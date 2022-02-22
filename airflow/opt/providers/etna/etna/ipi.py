from typing import List, Tuple, Dict

from airflow.utils.task_group import TaskGroup

from etna import (
    metis_etl,
    link,
    MetisEtlHelpers,
)
import re

from etna.etls.metis import MatchedRecordFolder
from etna.hooks.etna import File, Attribute

timepoint_regex = re.compile(r"^(?P<timepoint>MVIR1-HS\d+-DN?\d+)[A-Z]+.*")
patient_regex = re.compile(r"^(?P<patient>MVIR1-HS\d+)-D.*")
day_regex = re.compile(r".*-DN?(?P<day>\d+)$")
bulk_rna_seq_record_regex = re.compile(r"^bulk_RNASeq/raw/[^/]*")

pool_identifier_tri_regex = re.compile(r"^([^.]*\.[^.]*\.)(.*)$")
non_pool_identifier_tri_regex = re.compile("^([^_.]*[_.][^_.]*[_.])([^_.]*)(.*)$")


def correct_ipi_record_name(name):
    if "POOL" in name:
        name = name.replace("_", ".")
        match = pool_identifier_tri_regex.match(name)
        if match:
            name = match.groups(1) + match.groups(2).lower()
    else:
        match = non_pool_identifier_tri_regex.match(name)
        if match:
            name = (
                match.groups(1).replace("_", ".")
                + match.groups(2).lower()
                + match.groups(3)
            )
    return name


assert (
    correct_ipi_record_name("IPIPOOL001_P1_scRNA_XVIPBMC")
    == "IPIPOOL001.P1.scrna.xvipbmc"
)
assert (
    correct_ipi_record_name("IPILUNG094.N1.scrna.CD45_enriched")
    == "IPILUNG094.N1.scrna.CD45_enriched"
)
assert (
    correct_ipi_record_name("IPIPOOL001_P1_scRNA_XVIPBMC")
    == "IPIPOOL001.P1.scrna.xvipbmc"
)


def sc_file_linkers(helpers: MetisEtlHelpers, listed_matches):
    helpers.link_matching_file(
        listed_matches, "cite_antibody_key", re.compile(r"^ADT_keep_features\.list$")
    )
    helpers.link_matching_file(
        listed_matches, "mux_index_key", re.compile(r"^IDX_map\.tsx$")
    )
    helpers.link_matching_file(
        listed_matches,
        "processed_robject_rdata",
        re.compile(r".*_scTransformed_processed\.RData$"),
    )
    helpers.link_matching_file(
        listed_matches, "processed_umap", re.compile(r".*_umap\.pdf$")
    )
    helpers.link_matching_file(
        listed_matches, "processing_pipeline_parameters", re.compile(r"^cutoffs\.yml$")
    )
    helpers.link_matching_file(
        listed_matches,
        "filtered_counts_h5",
        re.compile(r"^filtered_feature_bc_matrix\.h5$"),
    )
    helpers.link_matching_file(
        listed_matches, "raw_counts_h5", re.compile(r"^raw_feature_bc_matrix\.h5$")
    )
    helpers.link_matching_file(
        listed_matches, "tenx_aligned_bam", re.compile(r"^possorted_genome_bam\.bam$")
    )
    helpers.link_matching_file(
        listed_matches,
        "tenx_aligned_bam_index",
        re.compile(r"^possorted_genome_bam\.bam\.bai$"),
    )
    helpers.link_matching_file(
        listed_matches, "tenx_cloupe_file", re.compile(r"^cloupe\.cloupe$")
    )
    helpers.link_matching_file(
        listed_matches, "tenx_metrics_csv", re.compile(r"^metrics_summary\.csv$")
    )
    helpers.link_matching_file(
        listed_matches, "tenx_molecule_info_h5", re.compile(r"^molecule_info\.h5$")
    )
    helpers.link_matching_file(
        listed_matches, "tenx_web_summary", re.compile(r"^web_summary\.html$")
    )


chem_version_regex = re.compile(r".*/(v1_chemistry|v2_chemistry|v3_chemistry)/.*")
prime_regex = re.compile(r".*/(3prime_|5prime_).*")


def link_specimen_and_chemistry(matches: List[MatchedRecordFolder]):
    for match in matches:
        biospecimen = match.record_name.split(".")[-1]
        prime = None
        version = "3"

        m = chem_version_regex.match(match.folder_path)
        if m:
            if "2" in m.group(1):
                version = "2"
            elif "1" in m.group(1):
                version = "1"

        m = prime_regex.match(match.folder_path)
        if m:
            if "3" in m.group(1):
                prime = "3"
            else:
                prime = "5"

        if prime is None:
            continue

        chemistry = f"10X_{prime}prime_v{version}"
        yield match, match.as_update(
            dict(
                biospecimen=biospecimen,
                chemistry=chemistry,
            )
        )


class MagmaRnaSeq:
    attributes: Dict[str, Attribute]
    raw: dict

    def __init__(self, raw: Dict, attributes: Dict[str, Attribute]):
        self.attributes = attributes
        self.raw = raw

    ATTRIBUTES_TO_SKIP = {
        "raw_base_count",
        "filtered_base_count",
        "raw_read_count",
        "filtered_read_count",
        "filter_passing_bases",
        "aligned_bases",
    }

    CONTROL_ATTRIBUTES_TO_SKIP = {
        "compartment",
    }


single_cell_root = r"^single_cell_[^/]*"
# Pools are nested
pool_root = r".*/IPIPOOL[^/]+/IPIPOOL[^/]+"
# Non pooled are not
non_pool_root = r".*/IPI((?!POOL)[^/])+"


@metis_etl("ipi", "data", 11)
def ipi_data_metis_etl(helpers: MetisEtlHelpers):
    with TaskGroup("sc_rna_seq_pool_raw"):
        matches = helpers.find_record_folders(
            "sc_rna_seq_pool",
            re.compile(f"{single_cell_root}/raw/{pool_root}"),
            corrected_record_name=correct_ipi_record_name,
        )
        matches = helpers.filter_by_timur(matches)
        listed_matches = helpers.list_match_folders(matches)
        helpers.link_matching_files(
            listed_matches,
            "raw_fastq",
            file_regex=re.compile(r".*\.fastq\.gz$"),
            dry_run=False,
        )

    with TaskGroup("sc_rna_seq_pool_processed"):
        matches = helpers.find_record_folders(
            "sc_rna_seq_pool",
            re.compile(f"{single_cell_root}/processed/{pool_root}"),
            corrected_record_name=correct_ipi_record_name,
        )
        matches = helpers.filter_by_timur(matches)
        matches = link(dry_run=False)(link_specimen_and_chemistry)(matches)
        listed_matches = helpers.list_match_folders(matches)
        sc_file_linkers(helpers, listed_matches)

    with TaskGroup("sc_rna_seq_raw"):
        matches = helpers.find_record_folders(
            "sc_rna_seq",
            re.compile(f"{single_cell_root}/raw/{non_pool_root}"),
            corrected_record_name=correct_ipi_record_name,
        )
        matches = helpers.filter_by_timur(matches)
        listed_matches = helpers.list_match_folders(matches)
        helpers.link_matching_files(
            listed_matches,
            "raw_fastq",
            file_regex=re.compile(r".*\.fastq\.gz$"),
            dry_run=False,
        )

    with TaskGroup("sc_rna_seq_processed"):
        matches = helpers.find_record_folders(
            "sc_rna_seq",
            re.compile(f"{single_cell_root}/processed/{pool_root}"),
            corrected_record_name=correct_ipi_record_name,
        )
        matches = helpers.filter_by_timur(matches)
        matches = link(dry_run=False)(link_specimen_and_chemistry)(matches)
        listed_matches = helpers.list_match_folders(matches)
        sc_file_linkers(helpers, listed_matches)

    with TaskGroup("rna_seq_bulk"):
        matches = helpers.find_record_folders(
            "-", re.compile(r"^bulkRNASeq/processed/[^/]*/results")
        )
        listed_matches = helpers.list_match_folders(matches)

        @link(dry_run=True)
        def process_rna_seq_attributes(listed_matches):
            for match, files in listed_matches:
                for file in files:
                    if file.file_name == "rnaseq_table.tsv":
                        pass

        process_rna_seq_attributes(listed_matches)
