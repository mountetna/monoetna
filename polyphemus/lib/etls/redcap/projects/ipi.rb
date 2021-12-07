# Define the new classes dynamically

# In order to get consistent date-shifting, we'll include
#   this module.
module PatientTable
  def self.included(base)
    base.class_eval do
      # We override these methods via the included hook,
      #   otherwise just including the PatientTable module
      #   means they will be overwritten by the default
      #   method definitions.
      def patient_id(magma_temp_id)
        magma_temp_id.split(/-/)[1].sub(/.clin/,'')
      end

      def offset_id(magma_record_name)
        patient_id(magma_record_name)
      end

      def patch(magma_record_name, record)
        set_record_patient(magma_record_name, record)
      end

      def set_record_patient(magma_record_name, record)
        record[:patient] = patient_id(magma_record_name)
      end
    end
  end
end

define_model("Demographic").class_eval do
  include PatientTable
  
  def patch(id, record)
    set_record_patient(id, record)

    age_fields = [ 'age at specimen collection', 'age at diagnosis' ]
    if age_fields.any? { |f| f == record[:name].downcase }
      record[:value] = [ 89, record[:value].to_i ].min.to_s if record[:value]
    end
  end
end

define_model("Treatment").class_eval do
  include PatientTable
end

define_model("Evaluation").class_eval do
  include PatientTable
end

define_model("Diagnostic").class_eval do
  include PatientTable
end

def config
  {
    models: {
      diagnostic: {
        each: [ "record" ],
        scripts: [
          {
            attributes: {
              name: {
                redcap_field: "hpv",
                value: "label"
              },
              value: "hpv",
              notes: "hpv_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "msi",
                value: "label"
                },
                value: "msi",
                notes: "msi_notes"
            }
          },
          {
            attributes: {
            name: {
            redcap_field: "baseline_ldh",
            value: "label"
            },
            value: "baseline_ldh",
            notes: "ldh_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "ih_mlha1",
                value: "label"
              },
              value: "ih_mlha1",
              notes: "mmr_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "ih_pms2",
                value: "label"
              },
              value: "ih_pms2",
              notes: "mmr_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "ih_msh2",
                value: "label"
              },
              value: "ih_msh2",
              notes: "mmr_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "ih_msh6",
                value: "label"
              },
              value: "ih_msh6",
              notes: "mmr_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "brca1",
                value: "label"
              },
              value: "brca1",
              notes: "brca1_mut_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "brca2",
                value: "label"
              },
              value: "brca2",
              notes: "brca2_mut_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "ras",
                value: "label"
              },
              value: "ras",
              notes: "ras_mut_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "kras",
                value: "label"
              },
              value: "kras",
              notes: "kras_mut_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "nras",
                value: "label"
              },
              value: "nras",
              notes: "nras_mut_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "braf",
                value: "label"
              },
              value: "braf",
              notes: "braf_mut_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "egfr",
                value: "label"
              },
              value: "egfr",
              notes: "egfr_mut_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "alk",
                value: "label"
              },
              value: "alk",
              notes: "alk_mut_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "ros1",
                value: "label"
              },
              value: "ros1",
              notes: "ros1_mut_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "chek1",
                value: "label"
              },
              value: "chek1",
              notes: "chek1_mut_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "chek2",
                value: "label"
              },
              value: "chek2",
              notes: "chek2_mut_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "palb2",
                value: "label"
              },
              value: "palb2",
              notes: "palb2_mut_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "atm",
                value: "label"
              },
              value: "atm",
              notes: "atm_mut_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "vhl",
                value: "label"
              },
              value: "vhl",
              notes: "vhl_mut_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "mlh1",
                value: "label"
              },
              value: "mlh1",
              notes: "mlh1_mut_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "pms2",
                value: "label"
              },
              value: "pms2",
              notes: "pms2_mut_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "msh2",
                value: "label"
              },
              value: "msh2",
              notes: "msh2_mut_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "msh6",
                value: "label"
              },
              value: "msh6",
              notes: "msh6_mut_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "number_of_lung_mets",
                value: "label"
              },
              value: {
                redcap_field: "number_of_lung_mets",
                value: "value",
                exists: true
              },
              notes: "met_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "number_of_liver_mets",
                value: "label"
              },
              value: {
                redcap_field: "number_of_liver_mets",
                value: "value",
                exists: true
              },
              notes: "met_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "number_of_brain_mets",
                value: "label"
              },
              value: {
                redcap_field: "number_of_brain_mets",
                value: "value",
                exists: true
              },
              notes: "met_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "number_of_bone_mets",
                value: "label"
              },
              value: {
                redcap_field: "number_of_bone_mets",
                value: "value",
                exists: true
              },
              notes: "met_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "number_of_lymphatic_mets",
                value: "label"
              },
              value: {
                redcap_field: "number_of_lymphatic_mets",
                value: "value",
                exists: true
              },
              notes: "met_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "number_of_skin_mets",
                value: "label"
              },
              value: {
                redcap_field: "number_of_skin_mets",
                value: "value",
                exists: true
              },
              notes: "met_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "number_of_other_mets",
                value: "label"
              },
              value: {
                redcap_field: "number_of_other_mets",
                value: "value",
                exists: true
              },
              notes: "met_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "clinical_t_stage",
                value: "label"
              },
              value: "clinical_t_stage",
              notes: "clin_t_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "clinical_n_stage",
                value: "label"
              },
              value: "clinical_n_stage",
              notes: "clin_n_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "clinical_m_stage_3",
                value: "label"
              },
              value: "clinical_m_stage_3",
              notes: "clin_m_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "pathologic_t_stage",
                value: "label"
              },
              value: "pathologic_t_stage",
              notes: "path_t_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "pathologic_n_stage",
                value: "label"
              },
              value: "pathologic_n_stage",
              notes: "path_n_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "pathologic_m_stage",
                value: "label"
              },
              value: "pathologic_m_stage",
              notes: "path_m_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "lymph_nodes_sampled",
                value: "label"
              },
              value: "lymph_nodes_sampled",
              notes: "ln_sampled_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "positive_lymph_nodes",
                value: "label"
              },
              value: "positive_lymph_nodes",
              notes: "ln_sampled_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "grade",
                value: "label"
              },
              value: "grade",
              notes: "grade_diff_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "residual_tumor_classification",
                value: "label"
              },
              value: "residual_tumor_classification",
              notes: "residual_tumor_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "residual_tumor_location",
                value: "label"
              },
              value: "residual_tumor_location",
              notes: "residual_tumor_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "sub_histology",
                value: "label"
              },
              value: "sub_histology",
              notes: "sub_hist_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "histology",
                value: "label"
              },
              value: "histology",
              notes: "hist_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "adren_inva",
                value: "label"
              },
              value: "adren_inva",
              notes: "adren_inva_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "child_pugh_score",
                value: "label"
              },
              value: "child_pugh_score",
              notes: "child_pugh_score_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "clin_figo_notes",
                value: "label"
              },
              value: "clin_figo_notes",
              notes: "clin_figo_notes_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "clinical_figo_stage",
                value: "label"
              },
              value: "clinical_figo_stage",
              notes: "clinical_figo_stage_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "clinical_stage",
                value: "label"
              },
              value: "clinical_stage",
              notes: "clinical_stage_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "extra_nodal_ext",
                value: "label"
              },
              value: "extra_nodal_ext",
              notes: "extra_nodal_ext_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "fuhr_grade",
                value: "label"
              },
              value: "fuhr_grade",
              notes: "fuhr_grade_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "hepatic_fibrosis",
                value: "label"
              },
              value: "hepatic_fibrosis",
              notes: "hepatic_fibrosis_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "lvi",
                value: "label"
              },
              value: "lvi",
              notes: "lvi_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "margin_status",
                value: "label"
              },
              value: "margin_status",
              notes: "margin_status_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "met_at_pres",
                value: "label"
              },
              value: "met_at_pres",
              notes: "met_at_pres_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "mets",
                value: "label"
              },
              value: "mets",
              notes: "mets_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "other_mut",
                value: "label"
              },
              value: "other_mut",
              notes: "other_mut_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "p13k",
                value: "label"
              },
              value: "p13k",
              notes: "p13k_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "pathologic_stage",
                value: "label"
              },
              value: "pathologic_stage",
              notes: "pathologic_stage_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "pdl1_staining",
                value: "label"
              },
              value: "pdl1_staining",
              notes: "pdl1_staining_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "perineural_invasion",
                value: "label"
              },
              value: "perineural_invasion",
              notes: "perineural_invasion_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "surg_marg",
                value: "label"
              },
              value: "surg_marg",
              notes: "surg_marg_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "tumor_necrosis",
                value: "label"
              },
              value: "tumor_necrosis",
              notes: "tumor_necrosis_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "tumor_size",
                value: "label"
              },
              value: "tumor_size",
              notes: "tumor_size_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "tumor_throm",
                value: "label"
              },
              value: "tumor_throm",
              notes: "tumor_throm_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "primary_diagnosis",
                value: "label"
                },
              value: "primary_diagnosis",
              notes: "primary_diagnosis_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "lymph_vas_invasion",
                value: "label"
                },
              value: "lymph_vas_invasion",
              notes: "lymph_vas_invasion_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "worst_pat_invasion",
                value: "label"
                },
              value: "worst_pat_invasion",
              notes: "worst_pat_invasion_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "depth_invasion",
                value: "label"
                },
              value: "depth_invasion",
              notes: "depth_invasion_notes"
            }
          }
        ]
      },
      evaluation: {
        each: [ "record" ],
        scripts: [
          {
            attributes: {
              name: {
                value: "text",
                text: "Diagnosis"
              },
              status: {
                value: "text",
                text: "Local Recurrence"
              },
              date: {
                redcap_field: "local_recurrence_date",
                exists: true,
                value: "value"
              }
            }
          },
          {
            attributes: {
              name: {
                value: "text",
                text: "Diagnosis"
              },
              status: {
                value: "text",
                text: "Regional Recurrence"
              },
              date: {
                redcap_field: "reg_recurrence_date",
                exists: true,
                value: "value"
              }
            }
          },
          {
            attributes: {
              name: {
                value: "text",
                text: "Diagnosis"
              },
              status: {
                value: "text",
                text: "Metastasis"
              },
              date: {
                redcap_field: "met_date",
                exists: true,
                value: "value"
              }
            }
          },
          {
            attributes: {
              name: {
                value: "text",
                text: "Diagnosis"
              },
              status: {
                value: "text",
                text: "Initial Diagnosis"

              },
              date: {
                redcap_field: "init_diag_date",
                exists: true,
                value: "value"
              },
              notes: "init_diag_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "followup",
                value: "label"
              },
              status: "followup",
              date: "followup_date",
              notes: "followup_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "death",
                value: "label"
              },
              status: "death",
              date: "death_date"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "response_neoadjuvant",
                value: "label"
              },
              status: "response_neoadjuvant",
              date: "response_date_neoadjuvant",
              notes: "response_date_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "response",
                value: "label"
              },
              status: "response",
              date: "response_date",
              notes: "response_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "tissue_extraction_site",
                value: "label"
              },
              status: "tissue_extraction_site",
              date: "tissue_extraction_date",
              notes: "tissue_extraction_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "os_time",
                value: "label"
              },
              status: "os_time"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "os",
                value: "label"
              },
              status: "os"
            }
          }
        ]
      },
      treatment: {
        each: [ "record", "repeat" ],
        scripts: [
          {
            comment: "drug treatments",
            each: [ "record", "repeat" ],
            attributes: {
              therapy: "treatment_type",
              category: {
                redcap_field: "treat_category",
                equals: "Drug",
                value: "value"
              },
              subtype: "drug_treat_type",
              name: "drug_name",
              value: "cycles",
              value_type: {
                value: "text",
                text: "cycles"
              },
              start_date: "regimen_start_date",
              stop_date: "regimen_stop_date",
              notes: "res_notes"
            }
          },
          {
            comment: "radiation treatments",
            each: [ "record", "repeat" ],
            attributes: {
              therapy: "treatment_type",
              category: {
                redcap_field: "treat_category",
                equals: "Radiation",
                value: "value"
              },
              subtype: "rad_treat_type",
              name: "brachytherapy",
              value: "fraction",
              value_type: {
                value: "text",
                text: "fraction"
              },
              start_date: "regimen_start_date",
              stop_date: "regimen_stop_date",
              notes: "res_notes"
            }
          },
          {
            comment: "procedure treatments, surgery",
            each: [ "record", "repeat" ],
            attributes: {
              therapy: "treatment_type",
              category: {
                redcap_field: "treat_category",
                equals: "Procedure",
                value: "value"
              },
              subtype: {
                redcap_field: "proc_type",
                equals: "Surgery",
                value: "value"
              },
              name: "proc_description",
              start_date: "regimen_start_date",
              stop_date: "regimen_stop_date",
              notes: "res_notes"
            }
          },
          {
            comment: "procedure treatments, biopsy",
            each: [ "record", "repeat" ],
            attributes: {
              therapy: "treatment_type",
              category: {
                redcap_field: "treat_category",
                equals: "Procedure",
                value: "value"
              },
              subtype: {
                redcap_field: "proc_type",
                equals: "Biopsy",
                value: "value"
              },
              name: "bio_description",
              start_date: "regimen_start_date",
              stop_date: "regimen_stop_date",
              notes: "res_notes"
            }
          }
        ]
      },
      demographic: {
        each: [ "record" ],
        scripts: [
          {
            attributes: {
              name: {
                redcap_field: "gender",
                value: "label"
              },
              value: {
                redcap_field: "gender",
                value: "value",
                exists: true
              },
              notes: "gender_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "age_at_diagnosis",
                value: "label"
              },
              value: {
                redcap_field: "age_at_diagnosis",
                value: "value",
                exists: true
              },
              notes: "age_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "age_specimen",
                value: "label"
              },
              value: {
                redcap_field: "age_specimen",
                value: "value",
                exists: true
              },
              notes: "age_notes_specimen"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "bmi",
                value: "label"
              },
              value: {
                redcap_field: "bmi",
                value: "value",
                exists: true
              },
              notes: "bmi_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "smoker",
                value: "label"
              },
              value: {
                redcap_field: "smoker",
                value: "value",
                exists: true
              },
              notes: "smoking_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "alcohol_use",
                value: "label"
              },
              value: {
                redcap_field: "alcohol_use",
                value: "value",
                exists: true
              },
              notes: "alcohol_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "race",
                value: "label"
              },
              value: {
                redcap_field: "race",
                value: "value",
                exists: true
              },
              notes: "race_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "ethnicity",
                value: "label"
              },
              value: {
                redcap_field: "ethnicity",
                value: "value",
                exists: true
              },
              notes: "ethnicity_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "prior_disease",
                value: "label"
              },
              value: {
                redcap_field: "prior_disease",
                value: "value",
                exists: true
              },
              notes: "prior_disease_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "prior_disease_2",
                value: "label"
              },
              value: {
                redcap_field: "prior_disease_2",
                value: "value",
                exists: true
              },
              notes: "prior_disease_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "prior_disease_3",
                value: "label"
              },
              value: {
                redcap_field: "prior_disease_3",
                value: "value",
                exists: true
              },
              notes: "prior_disease_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "prior_disease_4",
                value: "label"
              },
              value: {
                redcap_field: "prior_disease_4",
                value: "value",
                exists: true
              },
              notes: "prior_disease_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "prior_disease_5",
                value: "label"
              },
              value: {
                redcap_field: "prior_disease_5",
                value: "value",
                exists: true
              },
              notes: "prior_disease_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "prior_disease_6",
                value: "label"
              },
              value: {
                redcap_field: "prior_disease_6",
                value: "value",
                exists: true
              },
              notes: "prior_disease_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "prior_disease_7",
                value: "label"
              },
              value: {
                redcap_field: "prior_disease_7",
                value: "value",
                exists: true
              },
              notes: "prior_disease_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "prior_disease_8",
                value: "label"
              },
              value: {
                redcap_field: "prior_disease_8",
                value: "value",
                exists: true
              },
              notes: "prior_disease_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "prior_disease_9",
                value: "label"
              },
              value: {
                redcap_field: "prior_disease_9",
                value: "value",
                exists: true
              },
              notes: "prior_disease_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "prior_disease_10",
                value: "label"
              },
              value: {
                redcap_field: "prior_disease_10",
                value: "value",
                exists: true
              },
              notes: "prior_disease_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "primary_tumor_site",
                value: "label"
              },
              value: {
                redcap_field: "primary_tumor_site",
                value: "value",
                exists: true
              },
              notes: "primary_site_notes"
            }
          },
          {
            attributes: {
              name: {
                redcap_field: "met_site",
                value: "label"
              },
              value: {
                redcap_field: "met_site",
                value: "value",
                exists: true
              },
              notes: "met_site_notes"
            }
          }
        ]
      }
    }
  }
end
