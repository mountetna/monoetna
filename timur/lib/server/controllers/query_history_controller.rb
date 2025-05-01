class QueryHistoryController < Timur::Controller
  def create
    require_params(:query, :comment)

    query = QueryHistory.create(
      project_name: @params[:project_name],
      user: @user.email,
      query: @params[:query],
      comment: @params[:comment]
    )

    success_json(query: query.to_hash)
  end

  def list
    success_json(
      queries: QueryHistory.where(
        project_name: @params[:project_name],
        user: @user.email
      ).all.map(&:to_hash)
    )
  end

  def remove
    query = QueryHistory.where(
      project_name: @params[:project_name],
      id: @params[:id]
    ).first

    raise Etna::Forbidden unless @user.email == query.user

    query.delete

    success_json(success: true)
  end
end
