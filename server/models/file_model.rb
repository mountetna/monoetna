module FileModel

  class File < Sequel::Model
  end

  class Upload < Sequel::Model

    one_to_one :file
  end
end