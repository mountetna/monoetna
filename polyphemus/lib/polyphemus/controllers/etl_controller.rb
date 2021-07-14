require_relative 'controller'

class EtlController < Polyphemus::Controller
  def list
    success_json([ {
      project_name: @params[:project_name],
      etl: "redcap",
      name: "Redcap Loader",
      config: {}
    } ])
  end
end
