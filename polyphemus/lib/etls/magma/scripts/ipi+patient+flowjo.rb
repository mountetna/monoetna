return unless load_record_flow_jo('flojo_file')

add_stain_name('stain 1', 'treg')
add_stain_name('stain 2', 'nktb')
add_stain_name('stain 3', 'sort')
add_stain_name('stain 4', 'dc')
add_stain_name('stain 5', 'innate')

####### Edge Cases #######
# If not needed, keep the short versions
# If needed, swap for and modify the commented example
### Edge-case Load
def clean_name(name)
  puts "yeah yea yeah"
  name.gsub(/\s?,\s?/,',')
      .gsub(/ki67/i,'Ki67')
      .gsub(/foxp3/i,'FoxP3')
      .gsub(/PD-1/,'PD1') if name
end

####### Process the data #######
process_all_stains
process_all_populations
process_all_samples

