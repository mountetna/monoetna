require_relative 'controller'

class GnomonController < Magma::Controller
  def get
    grammar = require_grammar

    return success_json(grammar.to_hash)
  end

  def set
    require_param(:config, :comment)
    old_grammar = Magma::Gnomon::Grammar.for_project(project_name)

    version_number = (old_grammar&.version_number || 0) + 1

    errors = Magma::Gnomon::Grammar.validate(JSON.parse(@params[:config].to_json))

    return failure(422, errors: errors) unless errors.empty?

    grammar = Magma::Gnomon::Grammar.create(
      project_name: project_name,
      config: @params[:config],
      comment: @params[:comment],
      version_number: version_number
    )

    return success_json(grammar.to_hash)
  end

  def increment
    require_params(:project_name, :rule_name, :identifier_root)
    grammar = require_grammar
    rule = require_rule(grammar)

    next_value = rule.next(@params[:identifier_root])

    return success(next_value)
  rescue Magma::Gnomon::UnincrementableRuleError => e
    raise Etna::BadRequest, "Rule \"#{rule_name}\" is not incrementable."
  rescue Magma::Gnomon::UnrecognizedIdentifierError => e
    raise Etna::BadRequest, "Identifier root \"#{@params[:identifier_root]}\" does not match the rule definition for \"#{rule_name}\"."
  end

  def decompose
    grammar = require_grammar

    result = grammar.decompose(@params[:identifier])

    raise Etna::BadRequest, "Could not decompose identifier #{@params[:identifier]} for #{project_name}" unless result

    success_json(result)
  end

  def list
    grammar = require_grammar
    rule = require_rule(grammar)

    search_term = Regexp.new(@params[:regex] || ".*")
    identifiers = Magma::Gnomon::Identifier.where(
      project_name: project_name,
      rule: rule_name
    ).where { identifier =~ search_term }.map do |id|
      id.to_hash(@user)
    end

    success_json(identifiers)
  end

  def rule
    grammar = require_grammar
    rule = require_rule(grammar)

    rule_tokenization = rule.tokens.map{|t| grammar.tokens[t]&.merge(name: t) || { name: 'n', label: t } }
    success_json(rule: rule_tokenization)
  end

  def generate
    grammar = require_grammar
    rule = require_rule(grammar)

    raise Etna::BadRequest, "Invalid identifier \"#{@params[:identifier]}\" for rule \"#{rule_name}\" for #{project_name}" unless rule.valid?(@params[:identifier])

    identifier = Magma::Gnomon::Identifier.create(
      project_name: project_name,
      rule: rule_name,
      author: @user.display_name,
      identifier: @params[:identifier],
      grammar: grammar
    )

    success_json(identifier.to_hash(@user))
  rescue Sequel::UniqueConstraintViolation => e
    raise Etna::BadRequest, "Identifier \"#{@params[:identifier]}\" for rule \"#{rule_name}\" already exists"
  end

  def revisions
   grammars = Magma::Gnomon::Grammar.where(
     project_name: project_name
   ).reverse(:version_number).all

    success_json(grammars.map(&:to_revision))
  end

  private

  def require_grammar
    grammar = Magma::Gnomon::Grammar.for_project(project_name)

    if !grammar
      raise Etna::BadRequest, "No grammar found for project #{project_name}."
    end

    grammar
  end

  def require_rule(grammar)
    rule = grammar.rules[rule_name]

    raise Etna::BadRequest, "Unknown rule name, \"#{rule_name}\"." if rule.nil?

    rule
  end

  def project_name
    @params[:project_name]
  end

  def rule_name
    @params[:rule_name]
  end
end
