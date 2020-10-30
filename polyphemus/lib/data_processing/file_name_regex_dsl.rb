require_relative 'regex_dsl'

module FileNameRegexDsl
  def patient_name_regex(project_name = self.project_name)
    /#{project_name.upcase}\w+\d+/
  end

  def sample_name_regex(project_name = self.project_name)
    regex_chain(patient_name_regex(project_name), /\w+\d+/)
  end

  def flow_stain_name_regex(project_name = self.project_name)
    regex_chain(sample_name_regex(project_name), /_flow_\w+/)
  end

  def stain_name_regex(stain_names = @stain_to_names.values)
    regex_options(stain_names)
  end

  def self.included(other)
    unless other.included_modules.include?(RegexDsl)
      other.include(RegexDsl)
    end
  end
end