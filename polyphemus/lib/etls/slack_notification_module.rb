require_relative "../helpers"

module Polyphemus::WithSlackNotification
  include WithLogger

  def slack(title, channel, msg)
    logger.info(msg)
    `/bin/post-to-slack.sh "#{title}" "#{channel}" "#{msg}" || true`
  end
end
