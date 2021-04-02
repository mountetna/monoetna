return unless load_record_flow_jo('flojo_file_processed')

add_stain_name('stain 1', 'treg')
add_stain_name('stain 2', 'nktb')
add_stain_name('stain 3', 'sort')
add_stain_name('stain 4', 'dc')
add_stain_name('stain 5', 'innate')

def clean_name(name)
  name.gsub(/\s?,\s?/,',')
      .gsub(/ki67/i,'Ki67')
      .gsub(/foxp3/i,'FoxP3')
      .gsub(/PD-1/,'PD1') if name
end

process_all_populations

