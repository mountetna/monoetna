"""
Module container helper functions for processing metis files associated with both raw and processed single cell data,
both in pools and non pools.  Notably, these functions are versioned to account for changes in the pipeline that may not
necessaril ybe backported to all data sets.

Most of these helpers are accessible on the main MetisEtlHelper class as methods.
"""
from typing import Mapping, Optional, List, Tuple

import re

from airflow.exceptions import AirflowException
from airflow.models.xcom_arg import XComArg
from airflow.utils.task_group import TaskGroup
from etna.etls.metis import MatchedRecordFolder

from etna import MetisEtlHelpers, link

timepoint_regex = re.compile(r"^(?P<timepoint>[A-Z0-9]+-[A-Z0-9]+-[A-Z][A-Z]?\d+)[A-Z0-9]+.*")
patient_regex = re.compile(r"^(?P<patient>[A-Z0-9]+-[A-Z0-9]+)-.*")
day_regex = re.compile(r".*-[A-Z][A-Z]?(?P<day>\d+)$")

# TODO:  Would be great if a naming service provided this structure for us.
def create_and_link_parents(matches: List[MatchedRecordFolder]):
    for match in matches:
        if match.model_name == "sc_rna_seq_pool":
            yield match, match.as_update(dict(tube_name=match.record_name))

        if match.model_name != "sc_rna_seq":
            raise AirflowException(f"Model #{match.model_name} is not supported by this helper method.")

        timepoint_id = timepoint_regex.match(match.record_name)
        if not timepoint_id:
            raise AirflowException(f"Could not determine timpoint for sc_rna_seq record name #{match.record_name}, consider narrowing your identifier prefix")
        timepoint_id = timepoint_id.group("timepoint")

        patient_id = patient_regex.match(match.record_name)
        if not patient_id:
            raise AirflowException(f"Could not determine patient for sc_rna_seq record name #{match.record_name}, consider narrowing your identifier prefix")
        patient_id = patient_id.group("patient")

        day = day_regex.match(timepoint_id)
        if not day:
            raise AirflowException(f"Could not determine day for sc_rna_seq record name #{match.record_name}, consider narrowing your identifier prefix")

        day = int(day.group("day"))

        if "N" in timepoint_id:
            day *= -1

        update.update_record("patient", patient_id, dict())
        update.update_record(
            "timepoint", timepoint_id, dict(patient=patient_id, day=day)
        )
        update = match.as_update(dict(timepoint=timepoint_id))
        yield match, update

def link_single_cell_attribute_files_v1(
        helpers: MetisEtlHelpers,
        model_name: str,
        identifier_prefix: str,
        single_cell_root_prefix="single_cell_",
        attribute_linker_overrides: Optional[Mapping[str, str]]=None,
        dry_run=True,
) -> Tuple[XComArg, XComArg]:
    sc_file_linkers = dict(
        cite_antibody_key=r"^ADT_keep_features\.list$",
        mux_index_key=r"^IDX_map\.tsx$",
        processed_robject_rdata=r".*_scTransformed_processed\.RData$",
        processed_umap=r".*_umap\.pdf$",
        processing_pipeline_parameters=r"^cutoffs\.yml$",
        filtered_counts_h5=r"^filtered_feature_bc_matrix\.h5$",
        raw_counts_h5=r"^raw_feature_bc_matrix\.h5$",
        tenx_aligned_bam=r"^raw_feature_bc_matrix\.h5$",
        tenx_aligned_bam_index=r"^possorted_genome_bam\.bam\.bai$",
        tenx_cloupe_file=r"^cloupe\.cloupe$",
        tenx_metrics_csv=r"^metrics_summary\.csv$",
        tenx_molecule_info_h5=r"^molecule_info\.h5$",
        tenx_web_summary=r"^web_summary\.html$",
    )

    if attribute_linker_overrides:
        sc_file_linkers.update(**attribute_linker_overrides)

    with TaskGroup(f"{model_name}_processed"):
        processed_matches = helpers.find_record_folders(
            model_name,
            re.compile(f"^{single_cell_root_prefix}[^/]*/processed/{identifier_prefix}[^/]+"),
        )
        processed_matches = link(dry_run=dry_run)(create_and_link_parents)(processed_matches)
        listed_matches = helpers.list_match_folders(processed_matches)
        for attr, matcher in sc_file_linkers.items():
            helpers.link_matching_file(listed_matches, attr, re.compile(matcher), dry_run=dry_run)

    with TaskGroup(f"{model_name}_raw"):
        raw_matches = helpers.find_record_folders(
            model_name,
            re.compile(f"^{single_cell_root_prefix}[^/]*/raw/{identifier_prefix}[^/]+"),
        )

        raw_matches = link(dry_run=dry_run)(create_and_link_parents)(raw_matches)
        # create new records if needed
        listed_matches = helpers.list_match_folders(raw_matches)
        # make this a raw_fastq_files collection
        helpers.link_matching_files(
            listed_matches,
            "raw_fastq_files",
            file_regex=re.compile(r".*\.fastq\.gz$"),
            dry_run=dry_run,
        )

    return [processed_matches, raw_matches]
