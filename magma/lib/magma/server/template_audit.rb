require_relative 'controller'
require_relative '../template_audit'

class TemplateAuditController < Magma::Controller
  def report
    success_json(Magma::TemplateAudit.new.report)
  end
end
