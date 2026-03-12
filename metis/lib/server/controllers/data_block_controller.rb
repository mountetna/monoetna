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

  def vacuum_datablocks
    project_name = @params[:project_name]
    is_backfilled = (project_name == 'backfilled')
    include_projects = @params[:include_projects] || []
    include_projects = [include_projects] unless include_projects.is_a?(Array)
    
    # Parse commit parameter - defaults to false (dry-run mode for safety)
    commit = @params[:commit].nil? ? false : (@params[:commit] == true || @params[:commit].to_s.downcase == 'true')
    
    service = Metis::VacuumService.new(
      project_name: project_name,
      commit: commit,
      include_projects: include_projects,
      user: @user
    )
    
    response = service.vacuum_datablocks
    
    if commit
      event_log(
        event: 'vacuum_datablocks',
        message: is_backfilled ? "vacuumed backfilled orphaned datablocks" : "vacuumed orphaned datablocks for project #{project_name}",
        payload: {
          vacuumed_count: response[:summary][:total_vacuumed],
          space_freed: response[:summary][:space_freed]
        },
        consolidate: false
      )
    end
    
    success_json(response)
  end
end
