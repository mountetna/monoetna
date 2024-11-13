class Vesta
  class Project < Sequel::Model
    plugin :timestamps, update_on_create: true
    plugin :update_or_create  # https://sequel.jeremyevans.net/rdoc-plugins/classes/Sequel/Plugins/UpdateOrCreate.html


    def to_hash
      super.to_hash.yield_self do |h|
        h[:start_date] = h[:start_date].iso8601
        h[:created_at] = h[:created_at].iso8601
        h[:updated_at] = h[:updated_at].iso8601
        h[:principal_investigators] = h[:principal_investigators].to_a
        h[:data_types] = h[:data_types].to_a
        h.delete :id
        h
      end
    end
  end
end
