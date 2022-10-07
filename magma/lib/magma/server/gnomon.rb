require_relative 'controller'

class GnomonController < Magma::Controller
  def get
    grammar = Magma::Grammar.for_project(@params[:project_name])

    if !grammar
      raise Etna::BadRequest, "No grammar found for project #{@params[:project_name]}."
    end

    return success_json(grammar.to_hash)
  end

  def set
    require_param(:config)
    old_grammar = Magma::Grammar.for_project(@params[:project_name])

    version_number = (old_grammar&.version_number || 0) + 1

    errors = Magma::Grammar.validate(@params[:config])

    return failure(422, errors: errors) unless errors.empty?

    grammar = Magma::Grammar.create(
      project_name: @params[:project_name],
      config: @params[:config],
      version_number: version_number
    )

    return success_json(grammar.to_hash)
  end

  def increment
    require_params(:project_name, :rule_name, :identifier_root)

    grammar = Magma::Grammar.for_project(@params[:project_name])

    return failure(422, errors: ["No grammar defined for that project"]) if grammar.nil?

    rule = grammar.rules[@params[:rule_name]]

    return failure(422, errors: [
      "That rule is not incrementable."
    ]) unless rule.incrementable?

    return failure(422, errors: [
      "Identifier root \"#{@params[:identifier_root]}\" does not match the rule definition."
    ]) unless rule.regex(with_increment: false) =~ @params[:identifier_root]

    newest_identifier = Magma::Identifier.where(
      rule: rule.name
    ).where { identifier =~ RegExp.new(@params[:identifier_root]) }.order(:created_at).last

    decomposition = rule.decomposition(newest_identifier)
    decomposition["#{rule.name}_counter"] += 1

    next_identifier = Magma::Identifier.create(
      rule: rule.name,
      project_name: project_name,
      author: @user.name,
      identifier: rule.compose(decomposition)
    )

    return success_json(next_identifier.to_hash)
  end

  private

  def project_name
    @params[:project_name]
  end
end
