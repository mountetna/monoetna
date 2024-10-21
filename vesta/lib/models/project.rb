class Vesta
  class Project < Sequel::Model
    plugin :timestamps, update_on_create: true
    plugin :update_or_create  # https://sequel.jeremyevans.net/rdoc-plugins/classes/Sequel/Plugins/UpdateOrCreate.html


    def as_json
      JSON.parse(to_hash.to_json, symbolize_names: true)
    end
  end
end
