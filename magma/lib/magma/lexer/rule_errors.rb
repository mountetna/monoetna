class Magma
  class RuleParseError < Exception
  end

  class RecursiveRuleError < RuleParseError
  end

  class UndefinedRuleError < RuleParseError
  end
end
