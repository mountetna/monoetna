require_relative '../lib/etna/user'

describe Etna::User do
  it 'returns basic user info' do
    u = Etna::User.new(
      email: 'janus@two-faces.org',
      first: 'Janus',
      last: 'Bifrons',
      perm: 'a:labors;e:olympics,argo;v:constellations'
    )
    expect(u.first).to eq('Janus')
    expect(u.last).to eq('Bifrons')
    expect(u.email).to eq('janus@two-faces.org')
  end

  context "permissions" do
    before(:each) do
      @overlord = Etna::User.new(
        email: 'janus@two-faces.org',
        first: 'Janus',
        last: 'Bifrons',
        perm: 'a:administration'
      )
      @admin = Etna::User.new(
        email: 'polyphemus@etna.org',
        first: 'Polyphemus',
        last: 'Cyclops',
        perm: 'A:labors'
      )
      @editor = Etna::User.new(
        email: 'daedalus@two-faces.org',
        first: 'Daedalus',
        last: '',
        perm: 'E:labors'
      )
      @viewer = Etna::User.new(
        email: 'deino@graeae.org',
        first: 'Deino',
        last: 'Phorcides',
        perm: 'v:labors'
      )
    end

    it 'checks if the user can edit a project' do
      expect(@overlord.can_edit?(:labors)).to be_truthy
      expect(@admin.can_edit?(:labors)).to be_truthy
      expect(@editor.can_edit?(:labors)).to be_truthy
      expect(@viewer.can_edit?(:labors)).to be_falsy
    end
    it 'checks if the user can view a project' do
      expect(@overlord.can_view?(:labors)).to be_truthy
      expect(@admin.can_view?(:labors)).to be_truthy
      expect(@editor.can_view?(:labors)).to be_truthy
      expect(@viewer.can_view?(:labors)).to be_truthy
    end
    it 'checks if the user is an admin on the project' do
      expect(@overlord.is_admin?(:labors)).to be_truthy
      expect(@admin.is_admin?(:labors)).to be_truthy
      expect(@editor.is_admin?(:labors)).to be_falsy
      expect(@viewer.is_admin?(:labors)).to be_falsy
    end
    it 'checks if the user can see restricted data' do
      expect(@overlord.can_see_restricted?(:labors)).to be_falsy
      expect(@admin.can_see_restricted?(:labors)).to be_truthy
      expect(@editor.can_see_restricted?(:labors)).to be_truthy
      expect(@viewer.can_see_restricted?(:labors)).to be_falsy
    end
  end

  it "gives global permission to an administrator" do
    admin = Etna::User.new(
      email: 'janus@two-faces.org',
      first: 'Janus',
      last: 'Bifrons',
      perm: 'a:administration'
    )
    expect(admin.is_admin?(:administration)).to be_truthy
    expect(admin.is_admin?(:labors)).to be_truthy
  end
end
