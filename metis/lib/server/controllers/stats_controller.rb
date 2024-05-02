class StatsController < Metis::Controller
  def file_count_by_project
    projects_to_include = @params[:projects]

    data_blocks = Metis::DataBlock.left_join(:files, data_block_id: :id)
    data_blocks = data_blocks.where(project_name: projects_to_include) unless projects_to_include.nil?

    file_count_by_project = data_blocks
      .distinct(:data_block_id)
      .select_map(:project_name)
      .tally

    success_json(file_count_by_project)
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

    success_json(byte_count_by_project)
  end
end
