require_relative 'controller'

class LogController < Polyphemus::Controller
  MAX_MESSAGE_SIZE = 10_000

  def write
    require_params(:user, :event, :message)
    entry = Polyphemus::Log.create(
      application: @hmac.id,
      project_name: @params[:project_name],
      user: @params[:user],
      event: @params[:event],
      message: @params[:message][0..MAX_MESSAGE_SIZE],
      payload: @params[:payload],
      created_at: DateTime.now,
      hidden: false
    )

    success_json(log: entry.id)
  end

  def read
    query = {
      application: @params[:application],
      project_name: @params[:project_name] == 'administration' ? @params[:project_names] : @params[:project_name],
      user: @params[:user],
      event: @params[:event],
      message: @params[:message] ? Regexp.new(@params[:message]) : nil,
      created_at: @params[:from] && @params[:to] ? DateTime.parse(@params[:from])..DateTime.parse(@params[:to]) : nil,
      hidden: false,
    }.compact

    logs = Polyphemus::Log.where(
      query
    ).limit(100).all

    success_json(logs: logs.map(&:to_report))
  end

  def payload
    log = Polyphemus::Log[@params[:log_id].to_i]

    raise Etna::NotFound, "No such log #{@params[:log_id]}" unless log

    raise Etna::Forbidden, "Log id does not match #{@params[:project_name]}" unless log.project_name == @params[:project_name]

    success_json(payload: log.payload)
  end

  def hide
    log = Polyphemus::Log[@params[:log_id].to_i]

    raise Etna::NotFound, "No such log #{@params[:log_id]}" unless log

    log.hidden = true
    log.save

    success_json(log: @params[:log_id])
  end
end
