class QueryHistory < Sequel::Model
  def to_hash
    {
      id: id,
      query: query,
      user: user,
      project_name: project_name,
      comment: comment,
      created_at: created_at.iso8601
    }
  end
end
