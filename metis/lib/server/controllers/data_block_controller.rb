class DataBlockController < Metis::Controller
  def check
    raise Etna::BadRequest, "Improper md5!" unless @params[:md5s].all? do |md5|
      md5 =~ /^[a-f0-9]{32}$/i
    end
    found = Metis::DataBlock.where(md5_hash: @params[:md5s]).select_map(:md5_hash)
    success_json(found: found, missing: @params[:md5s] - found)
  end
end
