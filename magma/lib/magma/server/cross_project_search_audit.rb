require_relative 'controller'
require_relative '../audit/project_auditor'

class CrossProjectSearchAuditController < Magma::Controller
  def report
    success_json(Magma::ProjectAuditor.new(@project_name).report)
  end
end
