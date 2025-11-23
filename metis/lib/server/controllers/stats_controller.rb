class StatsController < Metis::Controller
  def file_count_by_project
    projects_to_include = @params[:projects]

    data_blocks = Metis::DataBlock.left_join(:files, data_block_id: :id)
    data_blocks = data_blocks.where(project_name: projects_to_include) unless projects_to_include.nil?

    file_count_by_project = data_blocks
      .distinct(:data_block_id)
      .select_map(:project_name)
      .tally

    success_json(count_detached!(file_count_by_project))
  end

  def file_count
    success_json(
      @params[:project_name] =>
      Metis::File
      .where(project_name: @params[:project_name])
      .distinct(:data_block_id).count
    )
  end

  def byte_count_by_project
    projects_to_include = @params[:projects]

    data_blocks = Metis::DataBlock.left_join(:files, data_block_id: :id)
    data_blocks = data_blocks.where(project_name: projects_to_include) unless projects_to_include.nil?

    byte_count_by_project = data_blocks
      .distinct(:data_block_id)
      .select_map([:size, :project_name])
      .group_by(&:last)
      .map { |project_name,sizes| [ project_name, sizes.map(&:first).sum ]}
      .to_h

    success_json(count_detached!(byte_count_by_project))
  end

  def byte_count
    success_json(
      @params[:project_name] => Metis::DataBlock.where(
        id: Metis::File.where(
          project_name: @params[:project_name]
        ).distinct(:data_block_id).select(:data_block_id)
      ).select_map(:size).sum
    )
  end

  def ledger
    verbose = @params[:verbose] && (@params[:verbose] == true || @params[:verbose].to_s.downcase == 'true')
    
    if @params[:legacy] && (@params[:legacy] == true || @params[:legacy].to_s.downcase == 'true')
      # Legacy mode: calculate vacuum stats for legacy datablocks
      orphaned_datablocks = Metis::DataBlockLedger.find_orphaned_datablocks_legacy
      
      vacuum_stats = {
        datablocks_can_vacuum: orphaned_datablocks.length,
        space_can_clear: orphaned_datablocks.sum(&:size)
      }
      
      if verbose
        vacuum_stats[:details] = Metis::DataBlockLedger.build_vacuum_details(orphaned_datablocks)
      end
      
      success_json(vacuum_stats)
    elsif @params[:project_name]
      # Project mode: high-level count summary + vacuum stats
      project_name = @params[:project_name]
      
      # Count all event types for this project
      event_counts = Metis::DataBlockLedger
        .where(project_name: project_name)
        .select_map(:event_type)
        .tally
      
      # Initialize all event types to 0 if they don't exist
      Metis::DataBlockLedger::EVENT_TYPES.each do |event_type|
        event_counts[event_type] ||= 0
      end
      
      # Calculate vacuum stats
      include_projects = @params[:include_projects] || []
      include_projects = [include_projects] unless include_projects.is_a?(Array)
      
      orphaned_datablocks = Metis::DataBlockLedger.find_orphaned_datablocks(project_name, include_projects: include_projects)
      
      vacuum_stats = {
        datablocks_can_vacuum: orphaned_datablocks.length,
        space_can_clear: orphaned_datablocks.sum(&:size),
        include_projects: include_projects
      }
      
      if verbose
        vacuum_stats[:details] = Metis::DataBlockLedger.build_vacuum_details(orphaned_datablocks, project_name)
      end
      
      success_json({
        project_name: project_name,
        event_counts: event_counts,
        vacuum: vacuum_stats
      })
    else
      raise Etna::BadRequest, "Must provide either 'project_name' or 'legacy' parameter"
    end
  end

  private

  def count_detached!(count_by_project)
    count_by_project["Detached"] = count_by_project.delete(nil)
    count_by_project.compact
  end
end
