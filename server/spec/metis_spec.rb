require './server/metis'

describe 'Integration' do
  
  before do

    @metis = Metis.new
  end

  describe Metis do

    it 'Implements the \'call\' method.' do

      expect(@metis).to respond_to(:call)
    end
  end
end