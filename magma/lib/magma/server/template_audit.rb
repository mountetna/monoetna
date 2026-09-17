require_relative 'controller'
require_relative '../template_audit'

class TemplateAuditController < Magma::Controller
  def models
    success_json(Magma::TemplateAudit.new.report)
  end
end
