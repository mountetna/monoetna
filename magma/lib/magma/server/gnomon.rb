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

  def decompose
    grammar = Magma::Gnomon::Grammar.for_project(@params[:project_name])

    result = grammar.decompose(@params[:identifier])

    raise Etna::BadRequest, "Could not decompose identifier #{@params[:identifier]} for #{@params[:project_name]}" unless result

    success_json(result)
  end
end
