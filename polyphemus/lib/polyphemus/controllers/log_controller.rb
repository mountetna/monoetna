require_relative 'controller'

class LogController < Polyphemus::Controller
  MAX_MESSAGE_SIZE = 10_000

  def write
    require_params(:user, :event, :message)

    if @params[:consolidate]
      entry = Polyphemus::Log.where(
        application: @hmac.id.to_s,
        project_name: @params[:project_name],
        user: @params[:user],
        event: @params[:event],
        message: @params[:message][0..MAX_MESSAGE_SIZE],
        hidden: false,
        created_at: (DateTime.now - 1)..(DateTime.now+1)
      ).first

      if entry
        entry.payload = merge_payload(
          entry.payload,
          JSON.parse(@params[:payload].to_json, symbolize_names: false)
        )
        entry.save

        return success_json(log: entry.id)
      end
    end

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

  private

  def merge_payload(p1, p2)
    p1.merge(p2) do |key, v1, v2|
      if v1.is_a?(Array) and v2.is_a?(Array)
        (v1 + v2)
      else
        v2
      end
    end
  end

  public

  def read
    query = {
      application: @params[:application],
      project_name: @params[:project_name] == 'administration' ? @params[:project_names] : @params[:project_name],
      user: @params[:user] ? Regexp.new(@params[:user]) : nil,
      event: @params[:event] ? Regexp.new(@params[:event]) : nil,
      message: @params[:message] ? Regexp.new(@params[:message]) : nil,
      created_at: @params[:from] && @params[:to] ? DateTime.parse(@params[:from])..DateTime.parse(@params[:to]) : nil,
      hidden: false,
    }.compact

    logs = Polyphemus::Log.where(
      query
    ).order(Sequel.desc(:created_at)).limit(100).all

    success_json(logs: logs.map(&:to_report))
  end

  def payload
    raise Etna::Forbidden unless @user.can_see_restricted?(@params[:project_name]) || @user.is_supereditor?

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
