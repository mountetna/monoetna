require_relative 'controller'

class FlagsController < Magma::Controller
  def get
    success_json({flags: @flags})
  end

  def set
    require_param(:flags)
    Magma::Flag.update_or_create(project_name, @params[:flags])
    success_json({success: 200})
  end

  def project_name
    @params[:project_name]
  end

end
