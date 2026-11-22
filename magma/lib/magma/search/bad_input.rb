class Magma
  module Search
    # Etna answers any error that is not an Etna::Error with a 500. Saving a
    # field or modality fails on ordinary bad input (blank or duplicate name,
    # unknown modality, bad kind), so report those as a 422 the editor can show.
    def self.bad_input_as_422
      yield
    rescue Sequel::ValidationFailed => e
      raise Etna::BadRequest, e.message.capitalize
    rescue Sequel::UniqueConstraintViolation
      raise Etna::BadRequest, 'That name already exists'
    rescue Sequel::ForeignKeyConstraintViolation
      raise Etna::BadRequest, 'That modality does not exist'
    rescue Sequel::CheckConstraintViolation
      raise Etna::BadRequest, "Type must be 'count' or 'filter'"
    end
  end
end
