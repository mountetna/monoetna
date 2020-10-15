require_relative '../helpers'
require_relative '../data_processing/xml_dsl'
require_relative '../data_processing/magma_dsl'
require_relative '../data_processing/file_name_regex_dsl'
require_relative '../data_processing/regex_dsl'
require_relative '../data_processing/project_metadata_dsl'

class EtlScriptRunner
  include XMLDsl
  include MagmaDsl
  include RegexDsl
  include FileNameRegexDsl
  include ProjectMetadataDsl

  attr_reader :script_name_parts
  def initialize(file_path)
    @file_path = File.expand_path(file_path)
    @script_name = File.basename(file_path)
    @script_name_parts = File.basename(@script_name, '.rb').split('+')
  end

  def project_name
    if @script_name_parts.length < 2
      raise "Improperly formatted script name #{@file_path}.  Script name should start with '<project name>_'."
    end

    @script_name_parts[0]
  end

  protected

  def run_script(binding = nil)
    File.open(@file_path, 'r') do |f|
      script = f.read
      eval(script, binding, @file_path)
    end
  end
end
