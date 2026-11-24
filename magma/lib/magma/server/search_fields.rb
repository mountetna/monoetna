require_relative 'controller'
require_relative '../search/bad_input'

class SearchFieldsController < Magma::Controller
  COLUMNS = %i[label query kind modality_id position].freeze

  def get
    success_json(fields: Magma::SearchField.eager(:modality).order(:position, :id).all.map(&:to_hash))
  end

  def create
    field = Magma::Search.bad_input_as_422 { Magma::SearchField.create(@params.slice(*COLUMNS)) }
    success_json(field: field.to_hash)
  end

  def set
    field = find_field
    Magma::Search.bad_input_as_422 { field.update(@params.slice(*COLUMNS)) }
    success_json(field: field.to_hash)
  end

  def delete
    find_field.destroy
    success_json(deleted: true)
  end

  private

  def find_field
    Magma::SearchField[@params[:id].to_i] || raise(Etna::NotFound, 'No such field')
  end
end
