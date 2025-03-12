require_relative 'controller'

class GnomonController < Magma::Controller
  def get
    grammar = require_grammar

    return success_json(grammar.to_hash)
  end

  def delete
    require_param(:identifiers)
    identifiers = Magma::Gnomon::Identifier.where(
      project_name: @params[:project_name],
      identifier: @params[:identifiers]
    ).select_map(:identifier)

    unknown = @params[:identifiers] - identifiers

    unless unknown.empty?
      raise Etna::BadRequest, "Unknown identifier(s): #{unknown.join(', ')}"
    end

    confirmation = Digest::MD5.hexdigest(@params[:identifiers].join)

    raise Etna::BadRequest, "Confirm deletion of #{@params[:identifiers].count} identifiers with code #{confirmation}" unless @params[:confirmation] == confirmation

    Magma::Gnomon::Identifier.where(
      project_name: @params[:project_name],
      identifier: @params[:identifiers]
    ).delete

    success_json(success: "Deleted #{@params[:identifiers].count} identifiers")
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
    ).where { identifier =~ search_term }.all

    if Magma.instance.get_project(project_name).models.has_key?(rule_name.to_sym)
      question = Magma::Question.new(
        project_name,
        [ rule_name, ["::identifier", "::in", identifiers.map(&:identifier)], "::all", "created_at" ],
        show_disconnected: false,
        restrict: !@user.can_see_restricted?(project_name),
        user: @user,
        timeout: Magma.instance.config(:query_timeout)
      )

      records_created_at = question.answer.to_h
    else
      records_created_at = {}
    end

    success_json(
      identifiers.map do |id|
        id.to_hash.merge(
          record_created_at: records_created_at[id.identifier]
        )
      end
    )
  end

  private
  def record(user)
    model = Magma.instance.get_model(project_name, rule)

    question = Magma::Question.new(
      project_name,
      record_query,
      show_disconnected: false,
      restrict: !user.can_see_restricted?(project_name),
      user: user,
      timeout: Magma.instance.config(:query_timeout)
    )

    question.answer
  end

  public
  def rule
    grammar = require_grammar
    rule = require_rule(grammar)

    rule_tokenization = rule.tokens.map do |t|
      grammar.parser.tokens.with_name(t) || { name: 'n', label: t }
    end
    success_json(rule: rule_tokenization)
  end

  def rules
    grammars = Magma::Gnomon::Grammar.where(project_name: @params[:project_names])
      .reverse(:project_name, :version_number).distinct(:project_name).all

    success_json(rules: grammars.to_h do |grammar|
      [
        grammar.project_name,
        grammar.parser.rules.to_h do |rule_name, rule|
          [ rule_name, rule.regex.source ]
        end
      ]
    end)
  end

  def project_rules
    grammar = Magma::Gnomon::Grammar.where(project_name: @params[:project_name])
      .reverse(:version_number).first

    success_json(rules: 
      grammar.parser.rules.to_h do |rule_name, rule|
        [ rule_name, rule.regex.source ]
      end
    )
  end

  def generate
    grammar = require_grammar
    rule = require_rule(grammar)

    raise Etna::BadRequest, "Invalid identifier \"#{@params[:identifier]}\" for rule \"#{rule_name}\" for #{project_name}" unless rule.valid?(@params[:identifier])

    decomposition = grammar.decompose(@params[:identifier])

    # find existing identifiers
    existing_identifiers = Magma::Gnomon::Identifier.where(
      project_name: project_name,
      identifier: decomposition[:rules].values.map{|r| r[:name]}
    ).all.map do |id|
      [ id.identifier, id ]
    end.to_h

    decomposition[:rules].each do |match_rule_name, rule|
      next if existing_identifiers[rule[:name]]

      rule_id = Magma::Gnomon::Identifier.create(
        project_name: project_name,
        rule: match_rule_name,
        author: @user.display_name,
        identifier: rule[:name],
        grammar: grammar
      )
    end

    created_identifiers = decomposition[:rules].map do |rule_name, rule|
      existing_identifiers[rule[:name]] ? nil : { rule_name: rule_name, name: rule[:name] }
    end.compact

    if created_identifiers
      Magma.instance.event_log(
        project_name: project_name,
        user: @user,
        event: 'create_name',
        message: "created #{created_identifiers.count} identifiers",
        payload: {
          identifiers: created_identifiers
        }
      )
    end

    success_json(decomposition)
  end

  def bulk_generate
    require_param(:names)
    grammar = require_grammar

    missing_rules = @params[:names].map {|n| n[:rule_name] }.uniq - grammar.parser.rules.rule_names

    unless missing_rules.empty?
      raise Etna::BadRequest, "Unknown rule names #{missing_rules.join(', ')}"
    end

    decompositions = @params[:names].map do |n|
      n.merge(
        decomposition: grammar.decompose(n[:name])
      )
    end

    bad_decompositions = decompositions.select do |n| n[:decomposition].nil? end

    unless bad_decompositions.empty?
      raise Etna::BadRequest, "Could not decompose names #{
        bad_decompositions.map{ |n| n[:name] }.join(', ')
      }"
    end

    identifier = nil

    decomposed_identifiers = decompositions.map do |d|
      d[:decomposition][:rules].values.map {|r| r[:name]}
    end.flatten.uniq

    # find existing identifiers
    existing_identifiers = Magma::Gnomon::Identifier.where(
      project_name: project_name,
      identifier: decomposed_identifiers
    ).all.map do |id|
      [ id.identifier, id ]
    end.to_h

    existing = []
    created = []

    created_identifiers = {}

    decompositions.each do |d|
      d[:decomposition][:rules].each do |match_rule_name, rule|
        next if created_identifiers[rule[:name]]

        if existing_identifiers[rule[:name]]
          existing.push(rule_name: match_rule_name, name: rule[:name] )
          next
        end

        rule_id = Magma::Gnomon::Identifier.create(
          project_name: project_name,
          rule: match_rule_name,
          author: @user.display_name,
          identifier: rule[:name],
          grammar: grammar
        )
        created.push(rule_name: match_rule_name, name: rule[:name] )

        created_identifiers[rule[:name]] = rule_id
      end
    end

    result = { created: created }
    result.update(existing: existing.uniq) unless existing.empty?

    unless created.empty?
      Magma.instance.event_log(
        project_name: project_name,
        user: @user,
        event: 'create_name',
        message: "created #{created.count} identifiers",
        payload: {
          identifiers: created
        }
      )
    end

    success_json(result)
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
    rule = grammar.parser.rules[rule_name]

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
