require_relative 'controller'

class GnomonController < Magma::Controller
  def get
    grammar = Magma::Gnomon::Grammar.for_project(@params[:project_name])

    if !grammar
      raise Etna::BadRequest, "No grammar found for project #{@params[:project_name]}."
    end

    return success_json(grammar.to_hash)
  end

  def set
    require_param(:config)
    old_grammar = Magma::Gnomon::Grammar.for_project(@params[:project_name])

    version_number = (old_grammar&.version_number || 0) + 1

    errors = Magma::Gnomon::Grammar.validate(JSON.parse(@params[:config].to_json))

    return failure(422, errors: errors) unless errors.empty?

    grammar = Magma::Gnomon::Grammar.create(
      project_name: @params[:project_name],
      config: @params[:config],
      version_number: version_number
    )

    return success_json(grammar.to_hash)
  end

  def increment
    require_params(:project_name, :rule_name, :identifier_root)

    grammar = Magma::Gnomon::Grammar.for_project(project_name)

    return failure(422, errors: ["No grammar defined for that project"]) if grammar.nil?

    rule = grammar.rules[@params[:rule_name]]

    return failure(422, errors: [
      "Unknown rule name, \"#{@params[:rule_name]}\"."
    ]) if rule.nil?

    next_value = rule.next(@params[:identifier_root])

    return success(next_value)
  rescue Magma::Gnomon::UnincrementableRuleError => e
    failure(422, errors: [
      "That rule is not incrementable."
    ])
  rescue Magma::Gnomon::UnrecognizedIdentifierError => e
    failure(422, errors: [
      "Identifier root \"#{@params[:identifier_root]}\" does not match the rule definition for \"#{@params[:rule_name]}\"."
    ])
  end

  def decompose
    grammar = Magma::Gnomon::Grammar.for_project(@params[:project_name])

    result = grammar.decompose(@params[:identifier])

    raise Etna::BadRequest, "Could not decompose identifier #{@params[:identifier]} for #{@params[:project_name]}" unless result

    success_json(result)
  end

  private

  def project_name
    @params[:project_name]
  end
end
