class Magma
  class RuleParser < RLTK::Parser

    def initialize(grammar)
      @grammar = grammar
      super()
    end

    production(:e) do
      @grammar['rules'].each do |token, rule|
        clause(rule) { |e| rule }
        clause(".e #{rule}") { |e| "#{@grammar['rule'][e.replace('.', '')]} #{rule}"}
      end
    end

    finalize
  end
end
