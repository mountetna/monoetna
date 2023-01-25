class DataBlockController < Metis::Controller
  def exists
    raise Etna::BadRequest, "Improper md5!" unless @params[:md5s].all? do |md5|
      md5 =~ /^[a-f0-9]{32}$/i
    end
    if @params[:project_name]
      found = Metis::DataBlock.where(md5_hash: @params[:md5s]).join(
        Metis::File.where(project_name: @params[:project_name]),
        data_block_id: :id
      ).select_map(:md5_hash)
    else
      found = Metis::DataBlock.where(md5_hash: @params[:md5s]).select_map(:md5_hash)
    end
    success_json(found: found, missing: @params[:md5s] - found)
  end
end
