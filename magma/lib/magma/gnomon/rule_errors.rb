class Magma
  module Gnomon
    class GrammarError < Exception
    end

    class RecursiveRuleError < GrammarError
    end

    class UndefinedRuleError < GrammarError
    end

    class UnrecognizedIdentifierError < GrammarError
    end

    class UnincrementableRuleError < GrammarError
    end
  end
end
