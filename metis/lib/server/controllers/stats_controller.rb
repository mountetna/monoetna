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
    if @params[:backfilled] && (@params[:backfilled] == true || @params[:backfilled].to_s.downcase == 'true')
      # Backfilled mode: calculate vacuum stats for backfilled datablocks
      orphaned_datablocks = Metis::DataBlockLedger.find_orphaned_datablocks_backfilled
      
      event_counts = Metis::DataBlockLedger.calculate_event_counts(nil)
      
      vacuum_stats = {
        datablocks_can_vacuum: orphaned_datablocks.length,
        space_can_clear: orphaned_datablocks.sum(&:size),
        details: Metis::DataBlockLedger.build_vacuum_details(orphaned_datablocks)
      }
      
      success_json({
        project_name: 'backfilled',
        event_counts: event_counts,
        vacuum: vacuum_stats
      })
    elsif @params[:project_name]
      # Tracked mode: high-level count summary + vacuum stats for tracked datablocks
      project_name = @params[:project_name]
      
      event_counts = Metis::DataBlockLedger.calculate_event_counts(project_name)
      
      # Calculate vacuum stats
      include_projects = @params[:include_projects] || []
      include_projects = [include_projects] unless include_projects.is_a?(Array)
      
      orphaned_datablocks = Metis::DataBlockLedger.find_orphaned_datablocks(project_name, include_projects: include_projects)
      
      project_breakdown = Metis::DataBlockLedger.calculate_project_breakdown(orphaned_datablocks, project_name, include_projects)
      
      vacuum_stats = {
        datablocks_can_vacuum: orphaned_datablocks.length,
        space_can_clear: orphaned_datablocks.sum(&:size),
        include_projects: include_projects,
        project_breakdown: project_breakdown,
        details: Metis::DataBlockLedger.build_vacuum_details(orphaned_datablocks, project_name)
      }
      
      success_json({
        project_name: project_name,
        event_counts: event_counts,
        vacuum: vacuum_stats
      })
    else
      raise Etna::BadRequest, "Must provide either 'project_name' (for tracked mode) or 'backfilled' parameter"
    end
  end

  private

  def count_detached!(count_by_project)
    count_by_project["Detached"] = count_by_project.delete(nil)
    count_by_project.compact
  end
end
