# Temporary
require 'securerandom'

class WorkflowError < Sequel::Model
  def self.mark_error!(hash:, message:)
    self.new(hash: hash, message: message, uuid: SecureRandom.uuid).tap(&:save)
  rescue Sequel::UniqueConstraintViolation
    find_error(hash: hash)
  end

  def self.find_error(hash:)
    self.find(hash: hash)
  end

  def clear!
    self.class.where(hash: self[:hash], uuid: self.uuid).delete
  end
end
