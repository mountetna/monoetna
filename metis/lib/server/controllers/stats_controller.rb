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
  private

  def count_detached!(count_by_project)
    count_by_project["Detached"] = count_by_project.delete(nil)
    count_by_project.compact
  end
end
