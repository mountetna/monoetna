class Vulcan
  class Snakemake
    class TargetParser
      def initialize(snakefile, config)
        @snakefile = snakefile # raw file
        @config = config  # yaml loaded
        @targets = {}  # Initialize the targets hash to store output targets as keys
      end

      # Public method to initiate parsing
      def parse
        preprocess
        @snakefile = filter_ui_rules
        parse_lines
        @targets  # Return the parsed targets with outputs as keys
      end

      private
      
      def preprocess
        # Use gsub! to replace only config variables whose values start with "output/"
        @snakefile.gsub!(/config\["([^"]+)"\]/) do
          key = Regexp.last_match(1)
          value = @config[key]
          # Only replace if the value starts with "output/"
          if value.nil?
            "config[\"#{key}\"]"  # Leave unchanged if the config key is not found
          elsif value.is_a?(String) && value.start_with?("output/")
            "\"#{value}\""  # Replace file paths starting with "output/"
          else
            "config[\"#{key}\"]"  # Leave unchanged for other config variables
          end
        end
      end

 
      # Removes rules that don't contain specified directives, like "ui rules"
      def filter_ui_rules
        directives = ["shell", "run", "script", "notebook", "wrapper", "cwl"]
        rules = @snakefile.split(/(?=rule\s+\w+:)/)
        filtered_rules = rules.select do |rule|
          directives.any? { |directive| rule.include?(directive + ":") }
        end
        filtered_rules.join("\n")
      end

      # Parses each line of the preprocessed file content
      def parse_lines
        current_rule = nil
        state = :none  # Possible states: :input, :output, :params, :none

        # Temporary variables for capturing rule details
        current_inputs = []
        current_outputs = []
        current_params = []

        @snakefile.each_line do |line|
          clean_line = clean_line(line)

          next if clean_line.empty?

          if rule_header?(clean_line)
            # When we encounter a new rule, store the previous rule (if any) and reset state
            store_target(current_outputs, current_inputs, current_params) if current_rule

            # Start a new rule
            current_rule = extract_rule_name(clean_line)
            current_inputs = []
            current_outputs = []
            current_params = []
            state = :none
            next
          end

          next unless current_rule

          case clean_line
          when /^input:/
            state = :input
          when /^output:/
            state = :output
          when /^params:/
            state = :params
          when /^(shell:|run:|rule )/
            state = :none
          else
            # Process line based on the current state
            process_line(clean_line, current_inputs, current_outputs, current_params, state)
          end
        end

        # Store the last parsed rule
        store_target(current_outputs, current_inputs, current_params) if current_rule
      end

      # Cleans a line by stripping whitespace and removing comments
      def clean_line(line)
        line.strip.gsub(/#.*$/, '')
      end

      # Checks if the line is a rule header
      def rule_header?(line)
        line.match(/^rule (\w+):/)
      end

      # Extracts the rule name from the rule header
      def extract_rule_name(line)
        line.match(/^rule (\w+):/)[1].to_sym
      end

      # Stores the parsed rule data in the @targets hash with output targets as keys
      def store_target(outputs, inputs, params)
        outputs.each do |output|
          @targets[output] = { inputs: inputs.uniq, params: params.uniq }
        end
      end

      # Processes a line based on the current state
      def process_line(line, inputs, outputs, params, state)
        case state
        when :input
          extract_io(line, inputs)
        when :output
          extract_io(line, outputs)
        when :params
          extract_params(line, params)
        end
      end

      # Extracts inputs or outputs from a line
      def extract_io(line, io_collection)
        # Handle key-value assignments (e.g., input1="file1", input2="file2")
        line.scan(/(\w+)\s*=\s*([^,]+)(?:,|$)/).each do |_, val|
          cleaned_val = clean_value(val)
          io_collection << cleaned_val unless io_collection.include?(cleaned_val)
        end

        # Handle direct list items (e.g., "file1", "file2")
        line.scan(/["']([^"']+)["']/).each do |match|
          io_collection << match.first unless io_collection.include?(match.first)
        end
      end

      # Extracts parameters from a line
      def extract_params(line, params)
        line.scan(/(\w+)\s*=/).map(&:first).each do |param|
          params << param unless params.include?(param)
        end
      end

      # Cleans a value by removing surrounding quotes if present
      def clean_value(val)
        val.strip!
        if val.start_with?('"') && val.end_with?('"') || val.start_with?("'") && val.end_with?("'")
          val[1..-2]
        else
          val
        end
      end
    end
  end
end
