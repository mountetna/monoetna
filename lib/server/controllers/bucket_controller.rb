class BucketController < Metis::Controller
  def list
    buckets = Metis::Bucket.where(
      project_name: @params[:project_name]
    ).all

    success_json(buckets: buckets.map(&:to_hash))
  end

  def create
  end

  def update
  end

  def remove
  end
end
