@flow_jo = load_flowjo(File.expand_path(__FILE__ + "/../patient-flojo_file-IPIBLAD001.wsp"))
@patient_record = load_magma_record("patient", "IPIBLAD001")

add_stain_name('stain 1', :treg)
add_stain_name('stain 2', :nktb)
add_stain_name('stain 3', :sort)
add_stain_name('stain 4', :dc)
add_stain_name('stain 5', :innate)

process_all_stains(@patient_record, @flow_jo)
pp @stain_panels

