class Vesta
  class Project < Sequel::Model
    plugin :timestamps, update_on_create: true
    plugin :update_or_create  # https://sequel.jeremyevans.net/rdoc-plugins/classes/Sequel/Plugins/UpdateOrCreate.html


  end
end
