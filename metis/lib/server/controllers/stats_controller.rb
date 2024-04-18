class StatsController < Metis::Controller
  def stats
    file_count_per_project = Metis::DataBlock
      .left_join(:files, data_block_id: :id)
      .distinct(:data_block_id)
      .select_map(:project_name)
      .tally

    byte_count_per_project = Metis::DataBlock
      .left_join(:files, data_block_id: :id)
      .distinct(:data_block_id)
      .select_map([:size, :project_name])
      .group_by(&:last)
      .map { |project_name,sizes| [ project_name, sizes.map(&:first).sum ]}
      .to_h

    stats_per_project = {}

    file_count_per_project.keys.each do |project|
      stats_per_project[project] = {
        file_count: file_count_per_project[project],
        byte_count: byte_count_per_project[project],
      }
    end

    success_json(stats_per_project)
  end
end
