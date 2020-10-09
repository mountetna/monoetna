return unless load_record_flow_jo('flojo_file')
add_stain_name('stain 1', 'treg')
add_stain_name('stain 2', 'nktb')
add_stain_name('stain 3', 'sort')
add_stain_name('stain 4', 'dc')
add_stain_name('stain 5', 'innate')

####### Edge Cases #######
# If not needed, keep the short versions
# If needed, swap for and modify the commented example

### Clean Load
# def clean_name name
#   name
# end

# def sample_name_from tube_name
#   if tube_name == Ipi::IPI.sample_name
#     return Regexp.last_match[0]
#   else
#     raise Magma::LoadFailed.new([ "Could not guess sample name from tube name '#{tube_name}'"])
#   end
# end

### Edge-case Load
def clean_name name
  name.gsub(/\s?,\s?/,',')
      .gsub(/ki67/i,'Ki67')
      .gsub(/foxp3/i,'FoxP3')
      .gsub(/PD-1/,'PD1') if name
end

def sample_name
  chain(:patient_name, /[A-Z][0-9]/)
end

def sample_name_from tube_name
  case tube_name
  when Ipi::IPI.sample_name
    return Regexp.last_match[0]
  when /[\W\_](?<code>[TN][0-9])/i
    return "#{@patient.ipi_number}.#{ Regexp.last_match[:code].upcase }"
  when /TUM(?:OR)?[\W\_]*(?<num>)[0-9]/i
    return "#{@patient.ipi_number}.T#{ Regexp.last_match[:num] }"
  when /TUM(?:OR)?[\W\_]*/i
    # just guess, least safe
    return "#{@patient.ipi_number}.T1"
  when /NORM(?:AL)?[\W\_]*(?<num>)[0-9]/i
    return "#{@patient.ipi_number}.N#{ Regexp.last_match[:num] }"
  when /NORM(?:AL)?[\W\_]*/i
    # just guess, least safe
    return "#{@patient.ipi_number}.N1"
  else
    raise Magma::LoadFailed.new([ "Could not guess sample name from tube name '#{tube_name}'"])
  end
end

####### Process the data #######

process_all_stains
process_all_populations

