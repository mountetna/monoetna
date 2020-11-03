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
      @ancestry ||= @parent ? [@parent.name] + @parent.ancestry : []
    end

    def statistics
      @statistics ||= @xml.css('> Subpopulations > Statistic').map do |stat|
        Statistic.new(stat)
      end
    end

    def to_hash
      {name => count}
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
        keyword_value(fluor_name.sub(/N$/, 'S')) if fluor_name
      end
    end

    def populations
      @populations ||= populations_for_node(@xml.css('SampleNode > Subpopulations > Population').first)
    end

    def populations_for_node node, parent = nil
      return [] unless node
      pop = Population.new(node, parent)
      children = node.css('> Subpopulations > Population').map do |child|
        populations_for_node(child, pop)
      end
      [pop, children].flatten
    end
  end

  class Group
    include XMLCensorInspect
    attr_reader :xml

    def initialize xml
      @xml = xml
    end

    def sample_refs
      @xml.css('SampleRefs > SampleRef').map { |sr| sr.attr("sampleID") }
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
      !sample.css(">SampleNode:named_like(\"sampleID\",\"^#{id}$\")", FlowJoXml::NameSearch.new).empty?
    end
    @samples[id] ||= FlowJoXml::Sample.new(match_samples.first)
  end

  def group name
    @groups[name] ||= FlowJoXml::Group.new(@xml.css("GroupNode:named_like(\"name\",\"#{name}\")", FlowJoXml::NameSearch.new))
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
      @flow_jo = FlowJoXml.new(file)
    end
  end

  def tubes_of_stain(flow_jo, stain)
    flow_jo.group(stain).sample_refs.map { |sr| flow_jo.sample(sr) }
  end

  def patient_flow_record_names(patient_record)
    query = [ 'flow',
              [ 'sample', 'patient', '::identifier', '::equals', patient_record['ipi_number'] ],
            '::all', '::identifier' ]

    request = Etna::Clients::Magma::QueryRequest.new(project_name: project_name, query: query)
    magma_client.query(request).answer.map { |flow| flow.last }.flatten
  end

  def process_all_populations(patient_record = @patient, flow_jo = @flow_jo)
    # Fetch all flow records for this patient, and pass that along to
    #   create_population_documents_for_stain.
    # If the flow record doesn't exist, throw an exception because the WSP
    #   is not well formed.

    existing_flow_record_names = patient_flow_record_names(patient_record)

    @stain_to_names.keys.each do |stain|
      create_population_documents_for_stain(flow_jo, stain, existing_flow_record_names)
    end
  end

  def flow_stain_name_from_tubename(tube_name)
    case tube_name
    when flow_stain_name_regex
      return tube_name.gsub('_', '.')
    else
      raise "Could not guess flow stain name from tube name '#{tube_name}', does not match #{flow_stain_name_regex.source}"
    end
  end

  def clean_name(name)
    name
  end

  def create_population_documents_for_stain(flow_jo, stain, existing_flow_record_names)
    tubes = tubes_of_stain(flow_jo, stain)

    if tubes.empty?
      return
    end

    tubes.each do |tube|
      flow_stain_name = flow_stain_name_from_tubename(tube.tube_name)

      puts "Processing flow #{flow_stain_name}."

      raise "WSP contains data for #{flow_stain_name}, but that flow record does not exist." unless existing_flow_record_names.include?(flow_stain_name)

      tube.populations.each do |pop|
        new_population(flow_stain_name, stain, pop)
      end
    end
  end

  def new_population(flow_stain_name, stain, pop)
    population_id = update_request.append_table("flow", flow_stain_name, "population", {
        ancestry: pop.ancestry.map { |name| clean_name(name) }.join("\t"),
        name: clean_name(pop.name),
        count: pop.count,
    })
  end
end