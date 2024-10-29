require 'json'

describe ExpungeController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  before(:each) do
  end

  def expunge(model_name, user_type=:editor, confirm_hash=nil)
    auth_header(user_type)
    json_post(:expunge, project_name: 'labors', model_name: model_name, confirm_hash: confirm_hash)
  end

  it 'fails for non-admins' do
    lion = create(:monster, name: 'Nemean Lion', species: 'lion')
    [:viewer, :guest, :editor].each do |role|
      expunge('monster', role)
      expect(last_response.status).to eq(403)
    end
  end

  it 'refuses to delete without a hash code' do
    lion = create(:monster, name: 'Nemean Lion', species: 'lion')
    hydra = create(:monster, name: 'Lernean Hydra', species: 'hydra')
    expunge('monster', :admin)
    expect(last_response.status).to eq(422)
    expect(Labors::Monster.count).to eq(2)
  end

  context 'with a hash code' do
    it 'deletes all records in a model' do
      lion = create(:monster, name: 'Nemean Lion', species: 'lion')
      hydra = create(:monster, name: 'Lernean Hydra', species: 'hydra')
      expunge('monster', :admin, Digest::MD5.hexdigest('labors.monster'))
      expect(last_response.status).to eq(200)
      expect(Labors::Monster.count).to eq(0)
      expect(json_body[:records_removed]).to eq(2)
    end

    it 'detaches, but does not delete, children of records from the expunged model' do
      labor = create(:labor, name: 'The Nemean Lion')
      lion = create(:monster, name: 'Nemean Lion', species: 'lion', labor: labor)
      expunge('labor', :admin, Digest::MD5.hexdigest('labors.labor'))
      expect(last_response.status).to eq(200)
      lion.refresh
      expect(Labors::Labor.count).to eq(0)
      expect(Labors::Monster.count).to eq(1)
      expect(lion.labor).to eq(nil)
    end
  end
end
