require_relative 'command'

# application.rb instantiates this for the project scoping.
# This generates a file, project-name.completion, which is sourced
# in build.sh to provide autocompletion in that environment.
module Etna
  class GenerateCompletionScript < Etna::Command
    def generate_for_command(command)
      completions = command.completions
      completions.each do |c|
        generate_start_match(c, false)
        write "fi"
        write "shift"
      end

      enable_flags(command.class)
      write 'while [[ "$#" != "0" ]]; do'
      generate_start_match([])
      generate_flag_handling

      write "else"
      write "return"
      write 'fi'
      write 'done'
      write "return"
    end

    def generate_flag_handling
      write %Q(elif [[ -z "$(echo $all_flag_completion_names | xargs)" ]]; then)
      write "return"
      write %Q(elif [[ "$all_flag_completion_names" =~ $1\\  ]]; then)
      write %Q(all_flag_completion_names="${all_flag_completion_names//$1\\ /}")
      write 'a=$1'
      write 'shift'
      write %Q(if [[ "$string_flag_completion_names" =~ $a\\  ]]; then)
      write 'if [[ "$#" == "1" ]];  then'
      write %Q(a="${a//--/}")
      write %Q(a="${a//-/_}")
      write %Q(i="_completions_for_$a")
      write %Q(all_completion_names="${!i}")
      write 'COMPREPLY=($(compgen -W "$all_completion_names" -- "$1"))'
      write 'return'
      write 'fi'
      write 'shift'
      write 'fi'
    end

    def generate_start_match(completions, include_flags=true)
      write 'if [[ "$#" == "1" ]];  then'
      write %Q(all_completion_names="#{completions.join(' ')}")
      write %Q(all_completion_names="$all_completion_names $all_flag_completion_names") if include_flags
      write %Q(if [[ -z "$(echo $all_completion_names | xargs)" ]]; then)
      write 'return'
      write 'fi'
      write 'COMPREPLY=($(compgen -W "$all_completion_names" -- "$1"))'
      write 'return'
    end

    def enable_flags(flags_container)
      boolean_flags = flags_container.boolean_flags
      string_flags = flags_container.string_flags
      flags = boolean_flags + string_flags
      write %Q(all_flag_completion_names="$all_flag_completion_names #{flags.join(' ')} ")
      write %Q(string_flag_completion_names="$string_flag_completion_names #{string_flags.join(' ')} ")

      string_flags.each do |flag|
        write %Q(declare _completions_for_#{flag_as_parameter(flag)}="#{completions_for(flag_as_parameter(flag)).join(' ')}")
      end
    end

    def generate_for_scope(scope)
      enable_flags(scope.class)
      write 'while [[ "$#" != "0" ]]; do'
      generate_start_match(scope.subcommands.keys)

      scope.subcommands.each do |name, command|
        write %Q(elif [[ "$1" == "#{name}" ]]; then)
        write 'shift'
        if command.class.included_modules.include?(CommandExecutor)
          generate_for_scope(command)
        else
          generate_for_command(command)
        end
      end

      generate_flag_handling

      write "else"
      write "return"
      write "fi"
      write 'done'
    end

    def program_name
      $PROGRAM_NAME
    end

    def execute
      name = File.basename(program_name)

      write <<-EOF
#!/usr/bin/env bash

function _#{name}_completions() {
  _#{name}_inner_completions "${COMP_WORDS[@]:1:COMP_CWORD}"
}

function _#{name}_inner_completions() {
  local all_flag_completion_names=''
  local string_flag_completion_names=''
  local all_completion_names=''
  local i=''
  local a=''
EOF
      generate_for_scope(parent)
      write <<-EOF
}

complete -o default -F _#{name}_completions #{name}
EOF

      File.open("#{name}.completion", 'w') { |f| f.write(@script) }
    end

    def write(string)
      @script ||= ""
      @script << string
      @script << "\n"
    end
  end
end
