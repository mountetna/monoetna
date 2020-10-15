require 'nokogiri'

module XMLDsl
  def load_xml(file_name)
    File.open(file_name) do |file|
      Nokogiri::XML(file.read)
    end
  end
end
