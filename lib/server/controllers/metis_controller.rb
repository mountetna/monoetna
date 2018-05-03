class Metis
  class Controller < Etna::Controller
    VIEW_PATH=File.expand_path('../views', __dir__)

    def success_json(obj)
      success(obj.to_json, 'application/json')
    end
  end
end
