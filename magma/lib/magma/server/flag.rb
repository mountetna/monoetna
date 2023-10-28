require_relative 'controller'

class FlagsController < Magma::Controller

  def initialize(request, action = nil)
    super
    @flags = Magma::Project.flags(@project_name)
  end

  def get
    success_json({flags: @flags})
  end

  def set
    require_param(:flags)
    @params[:flags].each do |flag|
      flag.each do |flag_name, flag_value|

        unless Magma::Flags.flag_valid?(flag_name.to_s, flag_value)
          raise Etna::BadRequest,  "Flag with name: \"#{flag_name}\", value: \"#{flag_value}\" is not registered in the Magma::Flags module."
        end

        update_key =  {project_name: @project_name, flag_name: flag_name.to_s}
        Magma::Flag.update_or_create(update_key, {value: flag_value})
      end
    end

    success_json({success: 200})
  end

end
