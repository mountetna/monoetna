require_relative 'controller'
require_relative '../search/bad_input'

class SearchModalitiesController < Magma::Controller
  COLUMNS = %i[name model position].freeze

  def get
    success_json(modalities: Magma::SearchModality.order(:position, :id).all.map(&:to_hash))
  end

  def create
    modality = Magma::Search.bad_input_as_422 { Magma::SearchModality.create(@params.slice(*COLUMNS)) }
    success_json(modality: modality.to_hash)
  end

  def set
    modality = find_modality
    Magma::Search.bad_input_as_422 { modality.update(@params.slice(*COLUMNS)) }
    success_json(modality: modality.to_hash)
  end

  def delete
    find_modality.destroy
    success_json(deleted: true)
  end

  private

  def find_modality
    Magma::SearchModality[@params[:id].to_i] || raise(Etna::NotFound, 'No such modality')
  end
end
