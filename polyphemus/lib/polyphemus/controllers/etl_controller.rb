require_relative 'controller'

# At base is the Job, a generic function provided as-is
#
# Each project may include any number of EtlConfig. Each
# EtlConfig configures a single type of Etl to run for
# this project.
#
# Each project may also configure any number of EtlRuns.
# Each EtlRun contains parameters describing the
# necessary runtime requirements (including secrets).
#
# Job
#   name
#   full_name
#   description
#   config_grammar
#   job_params
#
# JobConfig
#   project_name
#   job
#   config_name
#   config
#   run_every    # how often the run should be executed, -1 runs once, null disables
#   created_at
#   updated_at
#   ran_at       # when this job was last executed
#   completed_at # when this job last finished executing
#   output       # log of the output of this run
#   status       # the status of the job

class EtlController < Polyphemus::Controller
  def jobs
    return success_json(Polyphemus::Job.list.map(&:as_json))
  end

  def list
    success_json(
      Polyphemus::EtlConfig.exclude(archived: true).where(project_name: @params[:project_name]).all.map(&:as_json)
    )
  end

  def update
    etl_configs = Polyphemus::EtlConfig.exclude(archived: true).where(
      project_name: @params[:project_name],
      name: @params[:name]
    ).all

    if etl_configs.empty?
      raise Etna::FileNotFound, "No such etl #{@params[:name]} configured for project #{@params[:project_name]}"
    end

    if etl_configs.length > 1
      begin
        # repair broken configs
        *bad_etl_configs, etl_config = etl_configs.sort_by(&:updated_at)
        Polyphemus::EtlConfig.where(id: bad_etl_configs.map(&:id)).update(archived: true)
      end
    else
      etl_config = etl_configs.first
    end

    update = @params.slice(*(etl_config.columns - [:id, :project_name, :name]))

    if update[:config]
      # validate the configuration
      raise Etna::BadRequest, "Invalid configuration for etl \"#{etl_config.etl}\"" unless etl_config.validate_config(update[:config])
      new_etl_config = Polyphemus::EtlConfig.create(etl_config.as_json.merge(update))
      etl_config.update(archived: true, run_interval: Polyphemus::EtlConfig::RUN_NEVER)
      etl_config = new_etl_config
    else
      etl_config.update(update)
    end

    success_json(etl_config.as_json)
  end

  def create
    require_params(:job_type)

    etl_configs = Polyphemus::EtlConfig.exclude(archived: true).where(
      project_name: @params[:project_name],
      name: @params[:name]
    ).all

    unless etl_configs.empty?
      raise Etna::BadRequest, "There is already an etl #{@params[:name]} configured for project #{@params[:project_name]}"
    end

    unless Polyphemus::Job.from_name(@params[:job_type])
      raise Etna::BadRequest, "There is no such job type #{@params[:job_type]}"
    end

    etl_config = Polyphemus::EtlConfig.create(
      project_name: @params[:project_name],
      name: @params[:name],
      etl: @params[:job_type],
      config: {},
      run_interval: Polyphemus::EtlConfig::RUN_NEVER
    )

    success_json(etl_config.as_json)
  end

  def old_list
    success_json([ {
      project_name: @params[:project_name],
      etl: "redcap",
      name: "COMET Redcap Loader",
      ran_at: DateTime.parse("2021-03-04").iso8601,
      run: 86400,
      status: "completed",
      output: "This is the output log.\nIt has many lines",
      updated_at: DateTime.parse("2021-03-03").iso8601,
      created_at: DateTime.parse("2021-02-02").iso8601,
      config: {
        rna_seq: {
          each: [ "record", "event" ],
          invert: true,
          scripts: [
            {
              attributes: {
                tube_name: {
                  match: "-DN?\d+ETA",
                  value: "none",
                },
                biospecimen_date: "dos_esc_ta",
                biospecimen: {
                  value: "text",
                  text: "Endotracheal aspirate",
                },
              }
            },
            {
              attributes: {
                tube_name: {
                  match: "-DN?\d+PBMC",
                  value: "none",
                },
                biospecimen_date: "process_toc_4",
                biospecimen: {
                  value: "text",
                  text: "PBMCs",
                },
              }
            },
            {
              attributes: {
                tube_name: {
                  match: "-DN?\d+NAS",
                  value: "none",
                },
                biospecimen_date: "date_swab",
                biospecimen: {
                  value: "text",
                  text: "Nasal swab",
                },
              }
            },
          ],
        },
        sc_rna_seq: {
          each: [ "record", "event" ],
          invert: true,
          scripts: [
            {
              attributes: {
                tube_name: {
                  match: "-DN?\d+ETA",
                  value: "none",
                },
                biospecimen_date: "dos_esc_ta",
                biospecimen: {
                  value: "text",
                  text: "Endotracheal aspirate",
                },
              }
            },
            {
              attributes: {
                tube_name: {
                  match: "-DN?\d+PBMC",
                  value: "none",
                },
                biospecimen_date: "process_toc_4",
                biospecimen: {
                  value: "text",
                  text: "PBMCs",
                },
              }
            },
            {
              attributes: {
                tube_name: {
                  match: "-DN?\d+BLD",
                  value: "none",
                },
                biospecimen_date: "process_toc_4",
                biospecimen: {
                  value: "text",
                  text: "Whole blood",
                },
              }
            },
          ],
        },
        immunoassay: {
          each: [ "record", "event" ],
          invert: true,
          scripts: [
            {
              attributes: {
                name: {
                  match: "-DN?\d+PL",
                  value: "none",
                },
                biospecimen_date: "process_toc_4",
                biospecimen: {
                  value: "text",
                  text: "Plasma",
                },
              }
            },
            {
              attributes: {
                name: {
                  match: "-DN?\d+SR",
                  value: "none",
                },
                biospecimen_date: "process_toc_1",
                biospecimen: {
                  value: "text",
                  text: "Serum",
                },
              }
            },
          ],
        },
        cytof: {
          each: [ "record", "event" ],
          invert: true,
          scripts: [
            {
              attributes: {
                tube_name: {
                  match: "-DN?\d+BLD",
                  value: "none",
                },
                biospecimen_date: "process_toc_4",
                biospecimen: {
                  value: "text",
                  text: "Whole blood",
                },
              }
            },
          ],
        },
        patient: {
          each: [ "record", { "event" => "Enrollment" } ],
          scripts: [
            {
              attributes: {
                project: {
                  value: "text",
                  text: "COMET"
                },
                impacc_id: "impacc_id",
                comet_id: "id_study",
                d2b_id: "id_d2b",
                impacc_enrollment: "impacc_enroll",
                act1: "t15_act1_yn",
                act2: "t15_act2_yn",
                act3: "t15_act3_yn",
                consent: "sum_consent",
                admission_date: "hosp_dt",
                english_proficiency: "lep_yn",
                covid_status: "testing_status",
                hospital_site: "site_hospital",
                admission_level: "hosp_admitlevel",
                oxygen_homesupport: "hosp_homesupport",
                hosp_cause: "hosp_cause",
                hosp_cause_other: "hosp_cause_other",
                hosp_ed: "hosp_ed_yn",
                hosp_transfer: "hosp_transfer",
                race: "race",
                sex_at_birth: "gender",
                height: "height",
                weight: "weight",
                bmi: "vs_hosp_first_bmi",
                ethnicity: "ethnicity",
                age: "age_admit",
                deceased: "deceased",
                age_at_death: "age_death",
                discharge_date: "date_discharged",
                vent_duration: "vent_duration",
                standard_collection_date: "plan_date_standard",
                escalation_date: "plan_date_escalation",
                ta_collection_date: "plan_date_vent",
                hosp_los: "length_hosp",
                social_admit: "social_admit_yn",
                icu_los: "icu_los",
                admit_icu: "icu_admit_init",
                discharge_icu_to: "icu_disch_1",
                discharge_icu_date: "icu_disch_1_dt",
                hospitalization_outcome: "dispo",
                hospitalization_outcome_other: "dispo_other",
                dead_after_discharge: "dispo_late_death_yn",
                cause_of_death: "dispo_death_cause",
                cause_of_death_other: "dispo_death_cause_other",
                pulmonary_infection: "infect_pulm",
                non_pulmonary_infection: "infect_nonpulm",
                discharge_status: "discharge_status",
                admit_icu_date: "icu_admit_init_dt",
                readmit_icu: "icu_readmit",
                readmit_icu_date: "icu_readmit_dt_1",
                discharge_icu2_date: "icu_disch_2_dt",
                discharge_icu2_to: "icu_disch_2",
                readmit2_icu: "icu_readmit2",
                readmit_icu3_date: "icu_readmit_dt_2",
                discharge_icu3_date: "icu_disch_3_dt",
                discharge_icu3_to: "icu_disch_3",
                readmit3_icu: "icu_readmit3",
                readmit_icu4_date: "icu_readmit_dt_3",
                discharge_icu4_date: "icu_disch_4_dt",
                discharge_icu4_to: "icu_disch_4",
                primary_diag: "primary_diag",
                primary_diag_confirmed: "primary_diag_confirmed",
                primary_diag_subgroup: "primary_diag_subgroup",
                pulm_infection_timing: "infect_pulm_time",
                pneumonia_diagnosis_date: "infect_pulm_time_spec",
                bacteria_pneumonia_date: "infect_bact_time_spec",
                viral_pneumonia_date: "infect_viral_time_spec",
                viral_pneumonia_type: "viral_pneumonia_type",
                positive_culture: "infect_pulm_culture",
                pulm_pos_culture_source: "infect_pulm_culture_source",
                nonpulm_infection_timing: "infect_nonpulm_time",
                nonpulm_infection_date: "infect_nonpulm_time_spec",
                infect_nonpulm_micro: "infect_nonpulm_micro",
                nonpulm_pos_culture_source1: "infect_nonpulm_spec_source",
                infect_nonpulm_micro_2: "infect_nonpulm_micro2",
                nonpulm_pos_culture_source2: "infect_nonpulm_spec_src2",
                infect_nonpulm_micro_3: "infect_nonpulm_micro3",
                nonpulm_pos_culture_source3: "infect_nonpulm_spec_src3",
                infect_nonpulm_micro_4: "infect_nonpulm_micro4",
                nonpulm_pos_culture_source4: "infect_nonpulm_spec_src4",
                discharge_o2: "dispo_oxygen",
                oxygen_level: "dispo_oxygen_level",
                ards_berlin: "ards_berlin",
                ards_aecc: "ards_aecc",
                covid_pos: "covid_pos",
              }
            },
          ],
        },
        timepoint: {
          each: [ "record", "event" ],
          scripts: [
            {
              attributes: {
                o2therapy: "tx_o2therapy",
                o2flow: "tx_o2flow",
                o2current: "tx_o2_currentmode",
                o2high: "tx_o2mode",
                nippv: "tx_nippv",
                ecmo: "tx_ecmo",
                ecmo_start: "ecmo_start",
                ecmo_stop: "ecmo_stop",
                vasopressors_used: "tx_pressor",
                vasopressors_dose: "tx_pressor_dose",
                vasopressors: {
                  redcap_field: "tx_pressor_all",
                  value: "combine",
                  combine: ", "
                },
                rrt: "tx_rrt",
                rrt_type: "tx_rrt_type",
                sofa_resp: "sofa_resp",
                sofa_resp_calc: "sofa_resp_calc",
                sofa_coag: "sofa_coag",
                sofa_liver: "sofa_liver",
                sofa_cardio: "sofa_cardio",
                sofa_gcs: "sofa_gcs",
                sofa_renal: "sofa_renal",
                sofa_urine: "sofa_urine",
                sofa_score: "sofa_score",
                who_scale: "who_scale",
              }
            },
          ],
        },
        treatment: {
          each: [ "record", { "event" => "Enrollment" } ],
          scripts: [
            {
              attributes: {
                name: {
                  redcap_field: "meds_hosp_conplas",
                  equals: "Yes",
                  value: "label",
                },
                dose: "meds_hosp_conplas_dose",
                start: "meds_hosp_conplas_start",
                end: "meds_hosp_conplas_end",
                study: "meds_hosp_conplas_res",
                type: {
                  value: "text",
                  text: "COVID-19 targeted Medications",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "meds_hosp_hcq",
                  equals: "Yes",
                  value: "label",
                },
                dose: "meds_hosp_hcq_dose",
                start: "meds_hosp_hcq_start",
                end: "meds_hosp_hcq_end",
                study: "meds_hosp_hcq_res",
                type: {
                  value: "text",
                  text: "COVID-19 targeted Medications",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "meds_hosp_remdesivir",
                  equals: "Yes",
                  value: "label",
                },
                dose: "meds_hosp_remdesivir_dose",
                start: "meds_hosp_remdesivir_start",
                end: "meds_hosp_remdesivir_end",
                study: "meds_hosp_remdesivir_res",
                type: {
                  value: "text",
                  text: "COVID-19 targeted Medications",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "meds_hosp_toci",
                  equals: "Yes",
                  value: "label",
                },
                dose: "meds_hosp_toci_dose",
                start: "meds_hosp_toci_start",
                end: "meds_hosp_toci_end",
                study: "meds_hosp_toci_res",
                type: {
                  value: "text",
                  text: "COVID-19 targeted Medications",
                },
              }
            },
            {
              attributes: {
                name: {
                  value: "text",
                  text: "Chloroquine"
                },
                dose: {
                  redcap_field: "meds_hosp_chloro_dose",
                  value: "value",
                  exists: true
                },
                start: "meds_hosp_chloro_start",
                end: "meds_hosp_chloro_end",
                study: "meds_hosp_chloro_res",
                type: {
                  value: "text",
                  text: "COVID-19 targeted Medications",
                },
              }
            },
            {
              attributes: {
                name: {
                  value: "text",
                  text: "Lopinavir/ritonavir (Kaletra)"
                },
                dose: {
                  redcap_field: "meds_hosp_lopin_dose",
                  value: "value",
                  exists: true
                },
                start: "meds_hosp_lopin_start",
                end: "meds_hosp_lopin_end",
                study: "meds_hosp_lopin_res",
                type: {
                  value: "text",
                  text: "COVID-19 targeted Medications",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "meds_hosp_hiv_name",
                  exists: true,
                  value: "value",
                },
                dose: "meds_hosp_hiv_dose",
                start: "meds_hosp_hiv_start",
                end: "meds_hosp_hiv_end",
                study: "meds_hosp_hiv_res",
                type: {
                  value: "text",
                  text: "COVID-19 targeted Medications",
                },
              }
            },
            {
              attributes: {
                name: {
                  value: "text",
                  text: "Interferon alpha"
                },
                dose: {
                  redcap_field: "meds_hosp_inter_a_dose",
                  exists: true,
                  value: "value"
                },
                start: "meds_hosp_inter_a_start",
                end: "meds_hosp_inter_a_end",
                study: "meds_hosp_inter_a_res",
                type: {
                  value: "text",
                  text: "COVID-19 targeted Medications",
                },
              }
            },
            {
              attributes: {
                name: {
                  value: "text",
                  text: "Interferon beta"
                },
                dose: {
                  redcap_field: "meds_hosp_inter_b_dose",
                  exists: true,
                  value: "value"
                },
                start: "meds_hosp_inter_b_start",
                end: "meds_hosp_inter_b_end",
                study: "meds_hosp_inter_b_res",
                type: {
                  value: "text",
                  text: "COVID-19 targeted Medications",
                },
              }
            },
            {
              attributes: {
                name: {
                  value: "text",
                  text: "Ribavirin"
                },
                dose: {
                  redcap_field: "meds_hosp_riba_dose",
                  exists: true,
                  value: "value"
                },
                start: "meds_hosp_riba_start",
                end: "meds_hosp_riba_end",
                study: "meds_hosp_riba_res",
                type: {
                  value: "text",
                  text: "COVID-19 targeted Medications",
                },
              }
            },
            {
              attributes: {
                name: {
                  value: "text",
                  text: "Oseltamivir (Tamiflu)"
                },
                dose: {
                  redcap_field: "meds_hosp_tam_dose",
                  exists: true,
                  value: "value"
                },
                start: "meds_hosp_tam_start",
                end: "meds_hosp_tam_end",
                study: "meds_hosp_tam_res",
                type: {
                  value: "text",
                  text: "COVID-19 targeted Medications",
                },
              }
            },
            {
              attributes: {
                name: {
                  value: "text",
                  text: "Baloxavir"
                },
                dose: {
                  redcap_field: "meds_hosp_balo_dose",
                  exists: true,
                  value: "value"
                },
                start: "meds_hosp_balo_start",
                end: "meds_hosp_balo_end",
                study: "meds_hosp_balo_res",
                type: {
                  value: "text",
                  text: "COVID-19 targeted Medications",
                },
              }
            },
            {
              attributes: {
                name: {
                  value: "text",
                  text: "Sarulimab"
                },
                dose: {
                  redcap_field: "meds_hosp_saru_dose",
                  exists: true,
                  value: "value"
                },
                start: "meds_hosp_saru_start",
                end: "meds_hosp_saru_end",
                study: "meds_hosp_saru_res",
                type: {
                  value: "text",
                  text: "COVID-19 targeted Medications",
                },
              }
            },
            {
              attributes: {
                name: {
                  value: "text",
                  text: "Anakinra (Kinaret)"
                },
                dose: {
                  redcap_field: "meds_hosp_ana_dose",
                  exists: true,
                  value: "value"
                },
                start: "meds_hosp_ana_start",
                end: "meds_hosp_ana_end",
                study: "meds_hosp_ana_res",
                type: {
                  value: "text",
                  text: "COVID-19 targeted Medications",
                },
              }
            },
            {
              attributes: {
                name: {
                  value: "text",
                  text: "IV vitamin C"
                },
                dose: {
                  redcap_field: "meds_hosp_vitc_dose",
                  exists: true,
                  value: "value"
                },
                start: "meds_hosp_vitc_start",
                end: "meds_hosp_vitc_end",
                study: "meds_hosp_vitc_res",
                type: {
                  value: "text",
                  text: "COVID-19 targeted Medications",
                },
              }
            },
            {
              each: [ "record", { "event" => "Enrollment" }, { "field" => "meds_hosp_cs_iv_spec" } ],
              attributes: {
                name: {
                  redcap_field: "meds_hosp_cs_iv_spec",
                  equals: "IV methylprednisolone",
                  value: "value"
                },
                dose: "meds_hosp_cs_mp_dose",
                start: "meds_hosp_cs_mp_start",
                end: "meds_hosp_cs_mp_end",
                type: {
                  value: "text",
                  text: "IV Steroids",
                },
              }
            },
            {
              each: [ "record", { "event" => "Enrollment" }, { "field" => "meds_hosp_cs_iv_spec" } ],
              attributes: {
                name: {
                  redcap_field: "meds_hosp_cs_iv_spec",
                  equals: "IV dexamethasone",
                  value: "value"
                },
                dose: "meds_hosp_cs_dex_dose",
                start: "meds_hosp_cs_dex_start",
                end: "meds_hosp_cs_dex_end",
                type: {
                  value: "text",
                  text: "IV Steroids",
                },
              }
            },
            {
              each: [ "record", { "event" => "Enrollment" }, { "field" => "meds_hosp_cs_iv_spec" } ],
              attributes: {
                name: {
                  redcap_field: "meds_hosp_cs_iv_spec",
                  equals: "IV hydrocortisone",
                  value: "value",
                },
                dose: "meds_hosp_cs_hc_dose",
                start: "meds_hosp_cs_hc_start",
                end: "meds_hosp_cs_hc_end",
                type: {
                  value: "text",
                  text: "IV Steroids",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "meds_hosp_cs_inh",
                  equals: "Yes",
                  value: "label",
                },
                dose: "meds_hosp_cs_inh_dose",
                start: "meds_hosp_cs_inh_start",
                end: "meds_hosp_cs_inh_end",
                type: {
                  value: "text",
                  text: "Inhaled Steroids",
                },
              }
            },
            {
              each: [ "record", { "event" => "Enrollment" }, { "field" => "meds_hosp_cs_spec" } ],
              attributes: {
                name: "meds_hosp_cs_spec",
                dose: "",
                start: "",
                end: "",
                study: "",
                type: {
                  value: "text",
                  text: "Systemic Steroids",
                },
              }
            },
            {
              each: [ "record", { "event" => "Enrollment" }, { "field" => "meds_hosp_cs_po_spec" } ],
              attributes: {
                name: {
                  redcap_field: "meds_hosp_cs_po_spec",
                  equals: "Dexamethasone",
                  value: "value"
                },
                dose: "meds_hosp_cs_po_dex_dose",
                start: "meds_hosp_cs_po_dex_start",
                end: "meds_hosp_cs_po_dex_end",
                type: {
                  value: "text",
                  text: "Oral/Enteric Steroids",
                },
              }
            },
            {
              each: [ "record", { "event" => "Enrollment" }, { "field" => "meds_hosp_cs_po_spec" } ],
              attributes: {
                name: {
                  redcap_field: "meds_hosp_cs_po_spec",
                  equals: "Prednisone",
                  value: "value"
                },
                dose: "meds_hosp_cs_po_pred_dose",
                start: "meds_hosp_cs_po_pred_start",
                end: "meds_hosp_cs_po_pred_end",
                type: {
                  value: "text",
                  text: "Oral/Enteric Steroids",
                },
              }
            },
            {
              each: [ "record", { "event" => "Enrollment" }, { "field" => "meds_hosp_cs_po_spec" } ],
              attributes: {
                name: "meds_hosp_cs_po_oth",
                dose: "meds_hosp_cs_po_oth_dose",
                start: "meds_hosp_cs_po_oth_start",
                end: "meds_hosp_cs_po_oth_end",
                type: {
                  value: "text",
                  text: "Oral/Enteric Steroids",
                },
              }
            }
          ],
        },
        symptom: {
          each: [ "record", { "event" => "Enrollment" } ],
          scripts: [
            {
              attributes: {
                name: {
                  redcap_field: "symptom_fever",
                  value: "label",
                },
                present: {
                  redcap_field: "symptom_fever",
                  in: ["Yes", "No", "Unknown"],
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "symptom_chills",
                  value: "label",
                },
                present: {
                  redcap_field: "symptom_chills",
                  in: ["Yes", "No", "Unknown"],
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "symptom_shaking",
                  value: "label",
                },
                present: {
                  redcap_field: "symptom_shaking",
                  in: ["Yes", "No", "Unknown"],
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "symptom_cough",
                  value: "label",
                },
                present: {
                  redcap_field: "symptom_cough",
                  in: ["Yes", "No", "Unknown"],
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "symptom_cough_prod",
                  value: "label",
                },
                present: {
                  redcap_field: "symptom_cough_prod",
                  in: ["Yes", "No", "Unknown"],
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "symptom_haemoptysis",
                  value: "label",
                },
                present: {
                  redcap_field: "symptom_haemoptysis",
                  in: ["Yes", "No", "Unknown"],
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "symptom_sorethroat",
                  value: "label",
                },
                present: {
                  redcap_field: "symptom_sorethroat",
                  in: ["Yes", "No", "Unknown"],
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "symptom_rhinorrhea",
                  value: "label",
                },
                present: {
                  redcap_field: "symptom_rhinorrhea",
                  in: ["Yes", "No", "Unknown"],
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "symptom_wheeze",
                  value: "label",
                },
                present: {
                  redcap_field: "symptom_wheeze",
                  in: ["Yes", "No", "Unknown"],
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "symptom_chestpain",
                  value: "label",
                },
                present: {
                  redcap_field: "symptom_chestpain",
                  in: ["Yes", "No", "Unknown"],
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "symptom_myalgia",
                  value: "label",
                },
                present: {
                  redcap_field: "symptom_myalgia",
                  in: ["Yes", "No", "Unknown"],
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "symptom_arthralgia",
                  value: "label",
                },
                present: {
                  redcap_field: "symptom_arthralgia",
                  in: ["Yes", "No", "Unknown"],
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "symptom_fatigue",
                  value: "label",
                },
                present: {
                  redcap_field: "symptom_fatigue",
                  in: ["Yes", "No", "Unknown"],
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "symptom_dyspnea",
                  value: "label",
                },
                present: {
                  redcap_field: "symptom_dyspnea",
                  in: ["Yes", "No", "Unknown"],
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "symptom_edema",
                  value: "label",
                },
                present: {
                  redcap_field: "symptom_edema",
                  in: ["Yes", "No", "Unknown"],
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "symptom_nonamb",
                  value: "label",
                },
                present: {
                  redcap_field: "symptom_nonamb",
                  in: ["Yes", "No", "Unknown"],
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "symptom_headache",
                  value: "label",
                },
                present: {
                  redcap_field: "symptom_headache",
                  in: ["Yes", "No", "Unknown"],
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "symptom_confusion",
                  value: "label",
                },
                present: {
                  redcap_field: "symptom_confusion",
                  in: ["Yes", "No", "Unknown"],
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "symptom_seizure",
                  value: "label",
                },
                present: {
                  redcap_field: "symptom_seizure",
                  in: ["Yes", "No", "Unknown"],
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "symptom_syncope",
                  value: "label",
                },
                present: {
                  redcap_field: "symptom_syncope",
                  in: ["Yes", "No", "Unknown"],
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "symptom_anosmia",
                  value: "label",
                },
                present: {
                  redcap_field: "symptom_anosmia",
                  in: ["Yes", "No", "Unknown"],
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "symptom_abdpain",
                  value: "label",
                },
                present: {
                  redcap_field: "symptom_abdpain",
                  in: ["Yes", "No", "Unknown"],
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "symptom_nausea",
                  value: "label",
                },
                present: {
                  redcap_field: "symptom_nausea",
                  in: ["Yes", "No", "Unknown"],
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "symptom_diarrhea",
                  value: "label",
                },
                present: {
                  redcap_field: "symptom_diarrhea",
                  in: ["Yes", "No", "Unknown"],
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "symptom_conjunctivitis",
                  value: "label",
                },
                present: {
                  redcap_field: "symptom_conjunctivitis",
                  in: ["Yes", "No", "Unknown"],
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "symptom_bleeding",
                  value: "label",
                },
                present: {
                  redcap_field: "symptom_bleeding",
                  in: ["Yes", "No", "Unknown"],
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "symptom_earpain",
                  value: "label",
                },
                present: {
                  redcap_field: "symptom_earpain",
                  in: ["Yes", "No", "Unknown"],
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "symptom_chestwall",
                  value: "label",
                },
                present: {
                  redcap_field: "symptom_chestwall",
                  in: ["Yes", "No", "Unknown"],
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "symptom_skinrash",
                  value: "label",
                },
                present: {
                  redcap_field: "symptom_skinrash",
                  in: ["Yes", "No", "Unknown"],
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "symptom_skinulcer",
                  value: "label",
                },
                present: {
                  redcap_field: "symptom_skinulcer",
                  in: ["Yes", "No", "Unknown"],
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "symptom_lymphadenopathy",
                  value: "label",
                },
                present: {
                  redcap_field: "symptom_lymphadenopathy",
                  in: ["Yes", "No", "Unknown"],
                  value: "value",
                },
              }
            },
          ],
        },
        clinical_lab: {
          each: [ "record", "event" ],
          scripts: [
            {
              attributes: {
                name: {
                  redcap_field: "lab_ptt",
                  value: "label",
                },
                value: {
                  redcap_field: "lab_ptt",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "lab_ptt",
                  value: "note",
                },
                time: {
                  redcap_field: "da",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "lab_pt",
                  value: "label",
                },
                value: {
                  redcap_field: "lab_pt",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "lab_pt",
                  value: "note",
                },
                time: {
                  redcap_field: "da",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "lab_inr",
                  value: "label",
                },
                value: {
                  redcap_field: "lab_inr",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "lab_inr",
                  value: "note",
                },
                time: {
                  redcap_field: "da",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "lab_pct",
                  value: "label",
                },
                value: {
                  redcap_field: "lab_pct",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "lab_pct",
                  value: "note",
                },
                time: {
                  redcap_field: "da",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "lab_ldh",
                  value: "label",
                },
                value: {
                  redcap_field: "lab_ldh",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "lab_ldh",
                  value: "note",
                },
                time: {
                  redcap_field: "da",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "lab_ddimer",
                  value: "label",
                },
                value: {
                  redcap_field: "lab_ddimer",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "lab_ddimer",
                  value: "note",
                },
                time: {
                  redcap_field: "da",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "lab_ferritin",
                  value: "label",
                },
                value: {
                  redcap_field: "lab_ferritin",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "lab_ferritin",
                  value: "note",
                },
                time: {
                  redcap_field: "da",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "lab_platelets",
                  value: "label",
                },
                value: {
                  redcap_field: "lab_platelets",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "lab_platelets",
                  value: "note",
                },
                time: {
                  redcap_field: "da",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "lab_wbc",
                  value: "label",
                },
                value: {
                  redcap_field: "lab_wbc",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "lab_wbc",
                  value: "note",
                },
                time: {
                  redcap_field: "da",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "lab_trop",
                  value: "label",
                },
                value: {
                  redcap_field: "lab_trop",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "lab_trop",
                  value: "note",
                },
                time: {
                  redcap_field: "da",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "lab_na",
                  value: "label",
                },
                value: {
                  redcap_field: "lab_na",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "lab_na",
                  value: "note",
                },
                time: {
                  redcap_field: "da",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "lab_lactate",
                  value: "label",
                },
                value: {
                  redcap_field: "lab_lactate",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "lab_lactate",
                  value: "note",
                },
                time: {
                  redcap_field: "da",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "lab_k",
                  value: "label",
                },
                value: {
                  redcap_field: "lab_k",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "lab_k",
                  value: "note",
                },
                time: {
                  redcap_field: "da",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "lab_il6",
                  value: "label",
                },
                value: {
                  redcap_field: "lab_il6",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "lab_il6",
                  value: "note",
                },
                time: {
                  redcap_field: "da",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "lab_hgb",
                  value: "label",
                },
                value: {
                  redcap_field: "lab_hgb",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "lab_hgb",
                  value: "note",
                },
                time: {
                  redcap_field: "da",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "lab_hct",
                  value: "label",
                },
                value: {
                  redcap_field: "lab_hct",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "lab_hct",
                  value: "note",
                },
                time: {
                  redcap_field: "da",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "lab_hco3",
                  value: "label",
                },
                value: {
                  redcap_field: "lab_hco3",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "lab_hco3",
                  value: "note",
                },
                time: {
                  redcap_field: "da",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "lab_gluc",
                  value: "label",
                },
                value: {
                  redcap_field: "lab_gluc",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "lab_gluc",
                  value: "note",
                },
                time: {
                  redcap_field: "da",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "lab_esr",
                  value: "label",
                },
                value: {
                  redcap_field: "lab_esr",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "lab_esr",
                  value: "note",
                },
                time: {
                  redcap_field: "da",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "lab_crp",
                  value: "label",
                },
                value: {
                  redcap_field: "lab_crp",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "lab_crp",
                  value: "note",
                },
                time: {
                  redcap_field: "da",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "lab_cr",
                  value: "label",
                },
                value: {
                  redcap_field: "lab_cr",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "lab_cr",
                  value: "note",
                },
                time: {
                  redcap_field: "da",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "lab_ck",
                  value: "label",
                },
                value: {
                  redcap_field: "lab_ck",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "lab_ck",
                  value: "note",
                },
                time: {
                  redcap_field: "da",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "lab_bun",
                  value: "label",
                },
                value: {
                  redcap_field: "lab_bun",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "lab_bun",
                  value: "note",
                },
                time: {
                  redcap_field: "da",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "lab_bnp",
                  value: "label",
                },
                value: {
                  redcap_field: "lab_bnp",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "lab_bnp",
                  value: "note",
                },
                time: {
                  redcap_field: "da",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "lab_bili",
                  value: "label",
                },
                value: {
                  redcap_field: "lab_bili",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "lab_bili",
                  value: "note",
                },
                time: {
                  redcap_field: "da",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "lab_ast",
                  value: "label",
                },
                value: {
                  redcap_field: "lab_ast",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "lab_ast",
                  value: "note",
                },
                time: {
                  redcap_field: "da",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "lab_anc",
                  value: "label",
                },
                value: {
                  redcap_field: "lab_anc",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "lab_anc",
                  value: "note",
                },
                time: {
                  redcap_field: "da",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "lab_alt",
                  value: "label",
                },
                value: {
                  redcap_field: "lab_alt",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "lab_alt",
                  value: "note",
                },
                time: {
                  redcap_field: "da",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "lab_alb",
                  value: "label",
                },
                value: {
                  redcap_field: "lab_alb",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "lab_alb",
                  value: "note",
                },
                time: {
                  redcap_field: "da",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "lab_lymph",
                  value: "label",
                },
                value: {
                  redcap_field: "lab_lymph",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "lab_lymph",
                  value: "note",
                },
                time: {
                  redcap_field: "da",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "lab_wbc_dly",
                  value: "label",
                },
                value: {
                  redcap_field: "lab_wbc_dly",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "lab_wbc_dly",
                  value: "note",
                },
                time: {
                  redcap_field: "lab_wbc_diff_time",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "lab_hosp_img_dly",
                  value: "label",
                },
                value: {
                  redcap_field: "lab_hosp_img_dly",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "lab_hosp_img_dly",
                  value: "note",
                },
                time: {
                  redcap_field: "lab_hosp_img_dly_tm",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "lab_mac_hosp_dly",
                  value: "label",
                },
                value: {
                  redcap_field: "lab_mac_hosp_dly",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "lab_mac_hosp_dly",
                  value: "note",
                },
                time: {
                  redcap_field: "lab_mac_hosp_dly_tm",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "lab_bac_hosp_dly",
                  value: "label",
                },
                value: {
                  redcap_field: "lab_bac_hosp_dly",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "lab_bac_hosp_dly",
                  value: "note",
                },
                time: {
                  redcap_field: "lab_bac_hosp_dly_tm",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "lab_eo_hosp_dly",
                  value: "label",
                },
                value: {
                  redcap_field: "lab_eo_hosp_dly",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "lab_eo_hosp_dly",
                  value: "note",
                },
                time: {
                  redcap_field: "lab_eo_hosp_dly_tm",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "t15_lab_fibrin",
                  value: "label",
                },
                value: {
                  redcap_field: "t15_lab_fibrin",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "t15_lab_fibrin",
                  value: "note",
                },
                time: {
                  redcap_field: "da",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "vs_diastolic_high",
                  value: "label",
                },
                value: {
                  redcap_field: "vs_diastolic_high",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "vs_diastolic_high",
                  value: "note",
                },
                time: {
                  redcap_field: "dttm_vs_diastolic_high",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "vs_diastolic_low",
                  value: "label",
                },
                value: {
                  redcap_field: "vs_diastolic_low",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "vs_diastolic_low",
                  value: "note",
                },
                time: {
                  redcap_field: "dttm_vs_diastolic_low",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "vs_gcs",
                  value: "label",
                },
                value: {
                  redcap_field: "vs_gcs",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "vs_gcs",
                  value: "note",
                },
                time: {
                  redcap_field: "dttm_vs_gcs",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "vs_gcs_low",
                  value: "label",
                },
                value: {
                  redcap_field: "vs_gcs_low",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "vs_gcs_low",
                  value: "note",
                },
                time: {
                  redcap_field: "dttm_vs_gcs_low",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "vs_hr_high",
                  value: "label",
                },
                value: {
                  redcap_field: "vs_hr_high",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "vs_hr_high",
                  value: "note",
                },
                time: {
                  redcap_field: "dttm_vs_hr_high",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "vs_hr_low",
                  value: "label",
                },
                value: {
                  redcap_field: "vs_hr_low",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "vs_hr_low",
                  value: "note",
                },
                time: {
                  redcap_field: "dttm_vs_hr_low",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "vs_map",
                  value: "label",
                },
                value: {
                  redcap_field: "vs_map",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "vs_map",
                  value: "note",
                },
                time: {
                  redcap_field: "dttm_vs_map",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "vs_map_high",
                  value: "label",
                },
                value: {
                  redcap_field: "vs_map_high",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "vs_map_high",
                  value: "note",
                },
                time: {
                  redcap_field: "dttm_vs_map_high",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "vs_map_low",
                  value: "label",
                },
                value: {
                  redcap_field: "vs_map_low",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "vs_map_low",
                  value: "note",
                },
                time: {
                  redcap_field: "dttm_vs_map_low",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "vs_o2sat_high",
                  value: "label",
                },
                value: {
                  redcap_field: "vs_o2sat_high",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "vs_o2sat_high",
                  value: "note",
                },
                time: {
                  redcap_field: "dttm_vs_o2sat_high",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "vs_o2sat_low",
                  value: "label",
                },
                value: {
                  redcap_field: "vs_o2sat_low",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "vs_o2sat_low",
                  value: "note",
                },
                time: {
                  redcap_field: "dttm_vs_o2sat_low",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "vs_rr_spon_high",
                  value: "label",
                },
                value: {
                  redcap_field: "vs_rr_spon_high",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "vs_rr_spon_high",
                  value: "note",
                },
                time: {
                  redcap_field: "dttm_vs_rr_spon_high",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "vs_rr_spon_low",
                  value: "label",
                },
                value: {
                  redcap_field: "vs_rr_spon_low",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "vs_rr_spon_low",
                  value: "note",
                },
                time: {
                  redcap_field: "dttm_vs_rr_spon_low",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "vs_sbp_low",
                  value: "label",
                },
                value: {
                  redcap_field: "vs_sbp_low",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "vs_sbp_low",
                  value: "note",
                },
                time: {
                  redcap_field: "dttm_vs_sbp_low",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "vs_sbp_high",
                  value: "label",
                },
                value: {
                  redcap_field: "vs_sbp_high",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "vs_sbp_high",
                  value: "note",
                },
                time: {
                  redcap_field: "dttm_vs_sbp_high",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "vs_temp_high",
                  value: "label",
                },
                value: {
                  redcap_field: "vs_temp_high",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "vs_temp_high",
                  value: "note",
                },
                time: {
                  redcap_field: "dttm_vs_temp_high",
                  value: "value",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "vs_temp_low",
                  value: "label",
                },
                value: {
                  redcap_field: "vs_temp_low",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "vs_temp_low",
                  value: "note",
                },
                time: {
                  redcap_field: "dttm_vs_temp_low",
                  value: "value",
                },
              }
            },
          ],
        },
        comorbidity: {
          each: [ "record", { "event" => "Enrollment" } ],
          scripts: [
            {
              attributes: {
                name: {
                  redcap_field: "comorb_smoking",
                  value: "label",
                },
                value: "comorb_smoking",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_smoking_packs",
                  value: "label",
                },
                value: "comorb_smoking_packs",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_vaping",
                  value: "label",
                },
                value: "comorb_vaping",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_vaping_times",
                  value: "label",
                },
                value: "comorb_vaping_times",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_thc_cbd",
                  value: "label",
                },
                value: "comorb_thc_cbd",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_thc_cbd_times",
                  value: "label",
                },
                value: "comorb_thc_cbd_times",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_drugs",
                  value: "label",
                },
                value: "comorb_drugs",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_drugs_opiates",
                  value: "label",
                },
                value: "comorb_drugs_opiates",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_drugs_cocaine",
                  value: "label",
                },
                value: "comorb_drugs_cocaine",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_drugs_other",
                  value: "label",
                },
                value: "comorb_drugs_other",
              }
            },
            {
              attributes: {
                name: {
                  value: "text",
                  text: "Other illicit drugs used",
                },
                value: "comorb_drugs_other_spec",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_etoh",
                  value: "label",
                },
                value: "comorb_etoh",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_etoh_drinks",
                  value: "label",
                },
                value: "comorb_etoh_drinks",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_partum",
                  value: "label",
                },
                value: "comorb_partum",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_hiv",
                  value: "label",
                },
                value: "comorb_hiv",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_hiv_art_compliance",
                  value: "label",
                },
                value: "comorb_hiv_art_compliance",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_aids_cd4",
                  value: "label",
                },
                value: "comorb_aids_cd4",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_aids_viral_load",
                  value: "label",
                },
                value: "comorb_aids_viral_load",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_tb",
                  value: "label",
                },
                value: "comorb_tb",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_isaric_cardiac",
                  equals: "Yes",
                  value: "label",
                },
                value: "comorb_isaric_cardiac_spec",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_isaric_cardiac",
                  equals: "No",
                  value: "label",
                },
                value: {
                  value: "text",
                  text: "No",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_htn",
                  value: "label",
                },
                value: "comorb_htn",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_isaric_asthma",
                  value: "label",
                },
                value: "comorb_isaric_asthma",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_isaric_pulm",
                  value: "label",
                },
                value: "comorb_isaric_pulm",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_homerespsupp",
                  value: "label",
                },
                value: "comorb_homerespsupp",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_homeo2_amount",
                  value: "label",
                },
                value: "comorb_homeo2_amount",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_isaric_ckd",
                  value: "label",
                },
                value: "comorb_isaric_ckd",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_isaric_ckd_rrt",
                  value: "label",
                },
                value: "comorb_isaric_ckd_rrt",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_isaric_liver_mild",
                  value: "label",
                },
                value: "comorb_isaric_liver_mild",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_isaric_liver_mod",
                  value: "label",
                },
                value: "comorb_isaric_liver_mod",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_isaric_neuro",
                  equals: "Yes",
                  value: "label",
                },
                value: "comorb_isaric_neuro_spec",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_isaric_neuro",
                  equals: "No",
                  value: "label",
                },
                value: {
                  value: "text",
                  text: "No",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_isaric_neoplasm",
                  equals: "Yes",
                  value: "label",
                },
                value: "comorb_isaric_neo_spec",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_isaric_neoplasm",
                  equals: "No",
                  value: "label",
                },
                value: {
                  value: "text",
                  text: "No",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_isaric_heme",
                  equals: "Yes",
                  value: "label",
                },
                value: "comorb_isaric_heme_spec",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_isaric_heme",
                  equals: "No",
                  value: "label",
                },
                value: {
                  value: "text",
                  text: "No",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_isaric_obesity",
                  value: "label",
                },
                value: "comorb_isaric_obesity",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_isaric_dm_comp",
                  value: "label",
                },
                value: "comorb_isaric_dm_comp",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_isaric_dm_uncomp",
                  value: "label",
                },
                value: "comorb_isaric_dm_uncomp",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_isaric_rheum",
                  equals: "Yes",
                  value: "label",
                },
                value: "comorb_isaric_rheum_spec",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_isaric_rheum",
                  equals: "No",
                  value: "label",
                },
                value: {
                  value: "text",
                  text: "No",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_isaric_dementia",
                  value: "label",
                },
                value: "comorb_isaric_dementia",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_isaric_malnut",
                  value: "label",
                },
                value: "comorb_isaric_malnut",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_isaric_other",
                  equals: "Yes",
                  value: "label",
                },
                value: "comorb_isaric_other_spec",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_isaric_other",
                  equals: "No",
                  value: "label",
                },
                value: {
                  value: "text",
                  text: "No",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_bmt",
                  value: "label",
                },
                value: "comorb_bmt",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_transplant",
                  equals: "Yes",
                  value: "label",
                },
                value: "comorb_transplant_organ",
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "comorb_transplant",
                  equals: "No",
                  value: "label",
                },
                value: {
                  value: "text",
                  text: "No",
                },
              }
            },
          ],
        },
        admission_lab: {
          each: [ "record", { "event" => "Enrollment" } ],
          scripts: [
            {
              attributes: {
                name: {
                  redcap_field: "hosp_vs_hr",
                  value: "label",
                },
                value: {
                  redcap_field: "hosp_vs_hr",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "hosp_vs_hr",
                  value: "note",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "hosp_vs_rr_spon",
                  value: "label",
                },
                value: {
                  redcap_field: "hosp_vs_rr_spon",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "hosp_vs_rr_spon",
                  value: "note",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "hosp_vs_rr_vented",
                  value: "label",
                },
                value: {
                  redcap_field: "hosp_vs_rr_vented",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "hosp_vs_rr_vented",
                  value: "note",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "hosp_vs_sbp",
                  value: "label",
                },
                value: {
                  redcap_field: "hosp_vs_sbp",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "hosp_vs_sbp",
                  value: "note",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "hosp_vs_dpb",
                  value: "label",
                },
                value: {
                  redcap_field: "hosp_vs_dpb",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "hosp_vs_dpb",
                  value: "note",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "hosp_vs_o2sat",
                  value: "label",
                },
                value: {
                  redcap_field: "hosp_vs_o2sat",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "hosp_vs_o2sat",
                  value: "note",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "hosp_vs_map_calc",
                  value: "label",
                },
                value: {
                  redcap_field: "hosp_vs_map_calc",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "hosp_vs_map_calc",
                  value: "note",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "hosp_vs_map",
                  value: "label",
                },
                value: {
                  redcap_field: "hosp_vs_map",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "hosp_vs_map",
                  value: "note",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "hsp_vs_fio2_calc",
                  value: "label",
                },
                value: {
                  redcap_field: "hsp_vs_fio2_calc",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "hsp_vs_fio2_calc",
                  value: "note",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "o2amount",
                  value: "label",
                },
                value: {
                  redcap_field: "o2amount",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "o2amount",
                  value: "note",
                },
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: "hosp_temp",
                  value: "label",
                },
                value: {
                  redcap_field: "hosp_temp",
                  value: "value",
                  exists: true,
                },
                unit: {
                  redcap_field: "hosp_temp",
                  value: "note",
                },
              }
            },
          ],
        },
      }
    } ])
  end
end
