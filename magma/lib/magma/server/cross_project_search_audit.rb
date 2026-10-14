require_relative 'controller'
require_relative '../cross_project_search_auditor'

class CrossProjectSearchAuditController < Magma::Controller
  def report
    success_json(Magma::CrossProjectSearchAuditor.new.report)
  end
end
