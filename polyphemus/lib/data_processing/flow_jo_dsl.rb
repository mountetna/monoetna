require 'nokogiri'
require 'net/http'

module XMLCensorInspect
  def inspect
    "#<#{self.class} #{
    instance_variables.map do |var|
      next if var == :"@xml"
      "#{var}=#{instance_variable_get(var).inspect}"
    end.compact.join(", ")
    }>"
  end
end

class FlowJoXml
  include XMLCensorInspect

  class Statistic
    include XMLCensorInspect
    attr_reader :xml
    def initialize xml
      @xml = xml
    end

    def value
      @value ||= @xml.attr('value').to_f
    end

    def fluor
      @fluor ||= @xml.attr('id')
    end
  end

  class Population
    include XMLCensorInspect
    attr_reader :xml, :parent
    def initialize xml, parent
      @xml = xml
      @parent = parent
    end

    def name
      @name ||= @xml.attr('name')
    end

    def count
      @count ||= @xml.attr('count').to_i
    end

    def ancestry
      @ancestry ||= @parent ? [ @parent.name ] + @parent.ancestry : []
    end

    def statistics
      @statistics ||= @xml.css('> Subpopulations > Statistic').map do |stat|
        Statistic.new(stat)
      end
    end

    def to_hash
      { name => count }
    end
  end

  class Sample
    include XMLCensorInspect
    attr_reader :xml
    def initialize xml
      @xml = xml
      @stains = {}
    end

    def keyword_value name
      @xml.css("Keyword[name=\"#{name}\"]").map do |keyword|
        keyword.attr("value")
      end.first
    end

    def keyword_name value
      @xml.css("Keyword[value=\"#{value}\"]").map do |keyword|
        keyword.attr("name")
      end.first
    end

    def tube_name
      @tube_name ||= keyword_value("TUBE NAME")
    end

    def stain_for_fluor fluor
      @stains[fluor] ||= begin
        fluor_name = keyword_name(fluor)
        keyword_value(fluor_name.sub(/N$/,'S')) if fluor_name
      end
    end

    def populations
      @populations ||= populations_for_node(@xml.css('SampleNode > Subpopulations > Population').first)
    end

    def populations_for_node node, parent=nil
      return [] unless node
      pop = Population.new(node, parent)
      children = node.css('> Subpopulations > Population').map do |child|
        populations_for_node(child, pop)
      end
      [ pop, children ].flatten
    end
  end

  class Group
    include XMLCensorInspect
    attr_reader :xml
    def initialize xml
      @xml = xml
    end

    def sample_refs
      @xml.css('SampleRefs > SampleRef').map{|sr| sr.attr("sampleID") }
    end
  end

  class NameSearch
    def named_like(node_set, att, txt)
      node_set.find_all do |node|
        node[att] =~ /#{txt}/i
      end
    end
  end

  def initialize file
    @xml = Nokogiri::XML(file.read)
    @samples = {}
    @groups = {}
  end

  def sample id
    match_samples = @xml.css("Sample").select do |sample|
      !sample.css(">SampleNode:named_like(\"sampleID\",\"^#{id}$\")",FlowJoXml::NameSearch.new).empty?
    end
    @samples[id] ||= FlowJoXml::Sample.new(match_samples.first)
  end

  def group name
    @groups[name] ||= FlowJoXml::Group.new(@xml.css("GroupNode:named_like(\"name\",\"#{name}\")",FlowJoXml::NameSearch.new))
  end
end

module FlowJoDsl
  def add_stain_name(stain, name)
    @stain_to_names ||= {}
    @stain_to_names[stain] = name
  end

  def load_record_flow_jo(attribute_name)
    unless record.include?(attribute_name)
      raise "load_record_flow_jo cannot find #{attribute_name} on #{model_name}"
    end

    unless (value = record[attribute_name]).is_a?(Hash) && value.include?('url')
      logger.info("Record #{attribute_name} does not have a #{attribute_name} set, skipping flowjo processing.")
      return false
    end

    tmp = Tempfile.new
    begin
      uri = URI(value['url'])
      metis_client = Etna::Clients::Metis.new(host: uri.scheme + '://' + uri.host, token: magma_client.token)
      metis_client.download_file(value['url']) do |chunk|
        tmp.write chunk
      end
      load_flowjo(tmp.path)
    ensure
      tmp.close!
    end

    true
  end

  def load_flowjo(file_name)
    File.open(file_name) do |file|
      FlowJoXml.new(file)
    end
  end

  def new_stain_panel(patient_ipi_number, name, channels)
    @stain_panels ||= []
    @stain_panels <<  { patient_ipi_number: patient_ipi_number, name: name, channels: channels }
  end


  def process_all_stains(patient_record, flow_jo)
    @stain_to_names.keys.each do |stain|
      create_stain_document(patient_record, flow_jo, stain)
    end
  end

  def create_stain_document(patient_record, flow_jo, stain)
    tube = tube_of_stain(flow_jo, stain)
    name = @stain_to_names[stain]

    if !tube
      return
    end

    if !name
      raise "Name for stain #{stain} not configured, call add_stain_name('#{stain}', 'some-name')"
    end

    channels = []
    tube.keyword_value("$PAR").to_i.times do |n|
      # Each tube matches one stain. put a list together for each one.
      n = n + 1
      channels << {
          fluor: tube.keyword_value("$P#{n}N"),
          antibody: tube.keyword_value("$P#{n}S"),
          number: n,
      }
    end

    new_stain_panel(patient_record['ipi_number'], name, channels)
  end

  def tube_of_stain(flow_jo, stain)
    flow_jo.group(stain).sample_refs.map{|sr| flow_jo.sample(sr)}.first
  end
  
  def process_all_populations(patient_record, flow_jo)
    @stain_to_names.keys.each do |stain|
      create_population_document_per_stain(patient_record, flow_jo, stain)
    end
  end

  def create_population_document_per_stain patient_record, flow_jo, stain
    tube = tube_of_stain(flow_jo, stain)
    puts "Creating population #{tube} #{stain}"
    time = DateTime.now
    sample_name = sample_name_from(tube)
    populations = []
    mfis = []
    tube.populations.each do |pop|
      populations << {
        stain: @stain_to_names[stain].to_s,
        sample: sample_name,
        ancestry: pop.ancestry.map{|name| clean_name(name)}.join("\t"),
        name: clean_name(pop.name),
        count: pop.count}

      pop.statistics.each do |stat|
        mfis << {
          population: pop,
          name: clean_name(tube.stain_for_fluor(stat.fluor)),
          fluor: stat.fluor,
          value: stat.value}
      end
    end

    new_population_set(patient_record['ipi_number'], sample_name, populations, mfis)
  end

  def new_population_set(patient_ipi_number, name, channels)
    @population_set ||= []
    @population_set << {
      patient_ipi_number: patient_ipi_number,
      sample_name: sample_name,
      populations: populations,
      mfis: mfis }
  end

end