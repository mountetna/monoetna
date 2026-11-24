require_relative 'controller'
require_relative '../search/cross_project_search'

class CrossProjectSearchController < Magma::Controller
  def query
    success_json(search.query(@params[:queries]))
  end

  private

  def search
    Magma::CrossProjectSearch.new(@user)
  end
end
