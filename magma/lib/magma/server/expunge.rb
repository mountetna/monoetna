require_relative 'controller'

class ExpungeController < Magma::Controller
  def action
    require_param :model_name

    confirm_hash = Digest::MD5.hexdigest("#{@project_name}.#{@params[:model_name]}")
    raise Etna::BadRequest, "Please confirm with hash #{confirm_hash}" unless @params[:confirm_hash] == confirm_hash

    model = Magma.instance.get_model(@project_name, @params[:model_name])

    model.attributes.values.select(&:is_child_type?).each do |att|
      att.link_model.dataset.update(att.foreign_id => nil)
    end

    removed_count = model.dataset.delete
    
    return success_json(model_name: @params[:model_name], records_removed: removed_count)
  end
end
