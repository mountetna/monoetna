class Metis
  class Bucket < Sequel::Model
    one_to_many :files

    def allowed?(user)
      name == 'files'
    end
  end
end
