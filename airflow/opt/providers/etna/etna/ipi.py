from typing import List, Dict

from airflow.utils.task_group import TaskGroup

from etna import (
    metis_etl,
    link,
    MetisEtlHelpers, get_project_name,
)
import re

from etna.etls.metis import MatchedRecordFolder
from etna.hooks.etna import Attribute, UpdateRequest
import csv
import json
import os.path

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
            name = match.group(1) + match.group(2).lower()
    else:
        match = non_pool_identifier_tri_regex.match(name)
        if match:
            name = (
                    match.group(1).replace("_", ".")
                    + match.group(2).lower()
                    + match.group(3)
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

rna_seq_processed_file_linkers = dict(
    cram=r".*\.deduplicated\.cram$",
    cram_index=r".*\.deduplicated\.cram\.crai$",
    junction=r".*\.junction$",
    unmapped_fastqs=r".*\.fastq\.gz$",
)

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


class RnaSeqAttrTable:
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

    @property
    def tube_name(self):
        return RnaSeq.true_tube_name(self.raw[''])

    attr_name_mapping: Dict[str, str] = dict(
        cell_number="cell_count",
        duplication_pct="duplication_rate",
        ribosomal_read_count="reads",
        input_read_count="input_reads",
        uniq_map_read_count="uniq_map_reads",
        multimap_gt20_read_count="multimapp_gt20_reads",
        chimeric_read_count="chimeric_reads",
        chromosomal_read_count="chromosomal",
        mitochondrial_read_count="mitochondrial",
        median_5prime_bias="median_5p_bias",
        median_3prime_bias="median_3p_bias",
        eisenberg_score="EHK",
        expressed_eisenberg_genes="expressed_EHK_genes",
    )

    @property
    def raw_mean_length(self):
        parts = self.raw["raw_mean_length"].split(",")
        return (float(parts[0]) + float(parts[-1])) / 2

    @property
    def filtered_mean_length(self):
        parts = self.raw["filtered_mean_length"].split(",")
        return (float(parts[0]) + float(parts[-1])) / 2

    @property
    def compartment(self):
        value = self.raw['compartment']
        if value in self.valid_compartment_values:
            return value

        match = re.compile('(.*[^\d]+)(\d+)$').match(value)
        if match:
            maybe_value = match.group(1)
            if maybe_value in self.valid_compartment_values:
                return maybe_value

        return "other"

    @property
    def valid_compartment_values(self):
        return self.attributes['compartment'].validation['value']

    @property
    def as_revision(self):
        result = {}

        for attr_name in self.attributes.keys():
            if self.should_skip(attr_name):
                continue
            result[attr_name] = self.get_value_for_attr_name(attr_name)

        return result

    def should_skip(self, attr_name):
        if attr_name not in self.raw:
            return True

        if attr_name not in self.attributes:
            return True

        attr: Attribute = self.attributes[attr_name]
        if attr.hidden:
            return True

        if attr.restricted:
            return True

        if attr.read_only:
            return True

        if attr.link_model_name:
            return True

        if attr_name in self.ATTRIBUTES_TO_SKIP:
            return True

        if RnaSeq.is_control(self.tube_name) and attr_name in self.CONTROL_ATTRIBUTES_TO_SKIP:
            return True

        return False

    def get_value_for_attr_name(self, attr_name):
        orig_attr_name = attr_name
        if attr_name in self.attr_name_mapping:
            attr_name = self.attr_name_mapping[attr_name]

        if attr_name == "raw_mean_length":
            value = self.raw_mean_length
        elif attr_name == "filtered_mean_length":
            value = self.filtered_mean_length
        elif attr_name == "compartment":
            value = self.compartment
        elif attr_name == 'tube_name':
            value = self.tube_name
        else:
            value = self.raw[attr_name]

        if orig_attr_name in self.attributes:
            attr_type = self.attributes[orig_attr_name].attribute_type
            if attr_type == 'float':
                value = float(value)
            if attr_type == 'integer':
                value = int(value)

        return value


class RnaSeq:
    renames: Dict[str, str] = json.loads(open(
        os.path.join(os.path.dirname(__file__), 'ipi_bulk_rna_renames.json')
    ).read())

    @staticmethod
    def control_name_of_rna_seq_tube(tube_name: str) -> str:
        control, plate, *_ = tube_name.split(".")
        _, control_type, *_ = control.split("_")
        control_type: str
        plate: str

        if 'jurkat' in control_type.lower():
            return f"Control_Jurkat.{plate.capitalize()}"
        return f"Control_UHR.{plate.capitalize()}"

    @staticmethod
    def is_control(tube_name: str):
        return tube_name.lower().startswith('control_')

    @staticmethod
    def is_non_cancer_sample(tube_name: str):
        return 'NASH' in tube_name or 'NAFLD' in tube_name

    @staticmethod
    def true_tube_name(tube_name: str):
        if RnaSeq.is_control(tube_name):
            return RnaSeq.control_name_of_rna_seq_tube(tube_name)

        if tube_name in RnaSeq.renames:
            return RnaSeq.renames[tube_name]

        return tube_name


assert RnaSeq.true_tube_name('CONTROL_jurkat.plate1') == 'Control_Jurkat.Plate1'
assert RnaSeq.true_tube_name('CONTROL_uhr.plate2') == 'Control_UHR.Plate2'

single_cell_root = r"^single_cell_[^/]*"
# Pools are nested
pool_root = r".*/IPIPOOL[^/]+/IPIPOOL[^/]+"
# Non pooled are not
non_pool_root = r".*/IPI((?!POOL)[^/])+"


def process_gene_table(file_reader, gene_names, attr_name):
    update: UpdateRequest = UpdateRequest()
    for row in csv.DictReader(file_reader, delimiter='\t'):
        for key in row.keys():
            if key == 'gene_id':
                continue
            tube_name = RnaSeq.true_tube_name(key)

            if RnaSeq.is_non_cancer_sample(tube_name):
                continue

            revision = update.update_record('rna_seq', tube_name, {})
            matrix = revision.setdefault(attr_name, [0.0] * len(gene_names))
            if row['gene_id'] in gene_names and row[key] is not None:
                matrix[gene_names.index(row['gene_id'])] = float(row[key])

    return update


def ensure_update_empty(existing, revision, diff):
    if revision.get('identifier') in {'IPICRC041.T3.rna.live', 'IPILUNG008.T1.rna.live'}:
        return True
    return not diff

bulk_rna_seq_plate_regex = re.compile(r'(plate\d+)_.*')
def plate_name_of_bulk_rna_seq_folder(record_name: str) -> str:
    match = bulk_rna_seq_plate_regex.match(record_name)
    return match.group(1).capitalize()

assert plate_name_of_bulk_rna_seq_folder("plate123_rnaseq_new") == "Plate123"

@metis_etl("ipi", "data", 12)
def ipi_data_metis_etl(helpers: MetisEtlHelpers, tail_files):
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
        for attr, matcher in sc_file_linkers.items():
            helpers.link_matching_file(listed_matches, attr, re.compile(matcher), dry_run=False)

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
            re.compile(f"{single_cell_root}/processed/{non_pool_root}"),
            corrected_record_name=correct_ipi_record_name,
        )
        matches = helpers.filter_by_timur(matches)
        matches = link(dry_run=False)(link_specimen_and_chemistry)(matches)
        listed_matches = helpers.list_match_folders(matches)
        for attr, matcher in sc_file_linkers.items():
            helpers.link_matching_file(listed_matches, attr, re.compile(matcher), dry_run=False)

    with TaskGroup("rna_seq_bulk"):
        with TaskGroup("plate"):
            plate_matches = helpers.find_record_folders('rna_seq_plate', re.compile(r'^bulkRNASeq/processed/plate\d+[^/]*'), corrected_record_name=plate_name_of_bulk_rna_seq_folder)

            @link(dry_run=True, attribute_name='project', validate_record_update=ensure_update_empty)
            def create_plates_and_set_project_name(matches):
                for match in matches:
                    yield match, "UCSF Immunoprofiler"

            plate_matches = create_plates_and_set_project_name(plate_matches)

            with TaskGroup("non_control_rna_seq"):
                non_control_rna_seq_matches = helpers.find_record_folders('rna_seq', re.compile(r'^output/IPI[^/]*\.[A-Z]+\d\.[^/]*'),
                                                              corrected_record_name=RnaSeq.true_tube_name,
                                                              source=plate_matches)

                @link(dry_run=True, validate_record_update=ensure_update_empty)
                def link_sample_plate_and_rna_seq(matches):
                    for match in matches:
                        match: MatchedRecordFolder

                        if RnaSeq.is_non_cancer_sample(match.record_name):
                            continue

                        parts = match.record_name.split('.')
                        if len(parts) < 3:
                            continue
                        sample = '.'.join(parts[:2])
                        yield match, match.as_update(dict(
                            rna_seq_plate=match.match_parent.record_name,
                            sample=sample,
                        ))

                link_sample_plate_and_rna_seq(non_control_rna_seq_matches)
                listed_non_control = helpers.list_match_folders(non_control_rna_seq_matches)

                for attr, matcher in rna_seq_processed_file_linkers.items():
                    helpers.link_matching_file(listed_non_control, attr, re.compile(matcher), dry_run=True)

            with TaskGroup("control_rna_seq"):
                control_rna_seq_matches = helpers.find_record_folders('rna_seq', re.compile(r'^output/control_[^/]*', re.IGNORECASE),
                                                                          corrected_record_name=RnaSeq.control_name_of_rna_seq_tube,
                                                                          source=plate_matches)

                @link(dry_run=True, validate_record_update=ensure_update_empty)
                def link_plate_and_rna_seq(matches):
                    for match in matches:
                        yield match, match.as_update(dict(
                            rna_seq_plate=match.match_parent.record_name,
                        ))

                link_plate_and_rna_seq(control_rna_seq_matches)
                listed_control = helpers.list_match_folders(control_rna_seq_matches)

                for attr, matcher in rna_seq_processed_file_linkers.items():
                    helpers.link_matching_file(listed_control, attr,
                                               re.compile(matcher), dry_run=True)

        @link(dry_run=True, validate_record_update=ensure_update_empty)
        def process_gene_count_and_tpm_tables(files):
            file_regex = re.compile(r'^bulkRNASeq/processed/[^/]*/results/gene_(counts|tpm)_table\.tsv')
            with helpers.hook.magma() as magma:
                attributes = magma.retrieve(get_project_name(), 'rna_seq', hide_templates=False).models[
                    'rna_seq'].template.attributes
                gene_names = attributes['gene_counts'].validation['value']
            for metis_file in files:
                if not file_regex.match(metis_file.file_path):
                    continue
                attr_name: str = metis_file.file_name.replace("_table.tsv", "")

                with helpers.hook.metis() as metis:
                    with metis.open_file(metis_file) as file_reader:
                        yield metis_file, process_gene_table(file_reader, gene_names, attr_name)

        process_gene_count_and_tpm_tables(tail_files)

        @link(dry_run=True, validate_record_update=ensure_update_empty)
        def process_attr_table(files):
            file_regex = re.compile(r'^bulkRNASeq/processed/[^/]*/results/rnaseq_table\.tsv')
            with helpers.hook.magma() as magma:
                attributes = magma.retrieve(get_project_name(), 'rna_seq', hide_templates=False).models[
                    'rna_seq'].template.attributes
            with helpers.hook.metis() as metis:
                for metis_file in files:
                    if not file_regex.match(metis_file.file_path):
                        continue

                    update: UpdateRequest = UpdateRequest()
                    with metis.open_file(metis_file) as file_reader:
                        # Discard first line
                        file_reader.readline()

                        for row in csv.DictReader(file_reader, delimiter='\t'):
                            table = RnaSeqAttrTable(row, attributes)
                            if RnaSeq.is_non_cancer_sample(table.tube_name):
                                continue
                            update.update_record('rna_seq', table.tube_name, table.as_revision)
                    yield metis_file, update

        process_attr_table(tail_files)
