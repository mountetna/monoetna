require_relative '../lib/etna/user'

describe Etna::User do
  it 'returns basic user info with token' do
    u = Etna::User.new(
      {
        email: 'janus@two-faces.org',
        name: 'Janus Bifrons',
        perm: 'a:labors;e:olympics,argo;v:constellations'
      },
      'xyz123randomtoken'
    )
    expect(u.name).to eq('Janus Bifrons')
    expect(u.email).to eq('janus@two-faces.org')
    expect(u.token).to eq('xyz123randomtoken')
  end

  it 'returns basic user info without token param' do
    u = Etna::User.new(
      {
        email: 'janus@two-faces.org',
        name: 'Janus Bifrons',
        perm: 'a:labors;e:olympics,argo;v:constellations'
      }
    )
    expect(u.name).to eq('Janus Bifrons')
    expect(u.email).to eq('janus@two-faces.org')
    expect(u.token).to eq(nil)
  end

  context "permissions" do
    before(:each) do
      @superuser = Etna::User.new(
        email: 'janus@two-faces.org',
        name: 'Janus Bifrons',
        perm: 'a:administration'
      )
      @supereditor = Etna::User.new(
        email: 'portunus@two-faces.org',
        name: 'Portunus',
        perm: 'e:administration'
      )
      @superviewer = Etna::User.new(
        email: 'lar@two-faces.org',
        name: 'Lar Familiaris',
        perm: 'v:administration'
      )
      @admin = Etna::User.new(
        email: 'polyphemus@etna.org',
        name: 'Polyphemus',
        perm: 'A:labors'
      )
      @editor = Etna::User.new(
        email: 'daedalus@two-faces.org',
        name: 'Daedalus',
        perm: 'E:labors'
      )
      @viewer = Etna::User.new(
        email: 'deino@graeae.org',
        name: 'Deino Phorcides',
        perm: 'v:labors'
      )
      @guest = Etna::User.new(
        email: 'simon@earth.org',
        name: 'Simon the Greek',
        perm: 'g:labors'
      )
    end

    it 'checks if the user can edit a project' do
      expect(@superuser.can_edit?(:labors)).to be_truthy
      expect(@supereditor.can_edit?(:labors)).to be_truthy
      expect(@superviewer.can_edit?(:labors)).to be_falsy
      expect(@admin.can_edit?(:labors)).to be_truthy
      expect(@editor.can_edit?(:labors)).to be_truthy
      expect(@viewer.can_edit?(:labors)).to be_falsy
      expect(@guest.can_edit?(:labors)).to be_falsy
    end
    it 'checks if the user can view a project' do
      expect(@superuser.can_view?(:labors)).to be_truthy
      expect(@supereditor.can_view?(:labors)).to be_truthy
      expect(@superviewer.can_view?(:labors)).to be_truthy
      expect(@admin.can_view?(:labors)).to be_truthy
      expect(@editor.can_view?(:labors)).to be_truthy
      expect(@viewer.can_view?(:labors)).to be_truthy
      expect(@guest.can_view?(:labors)).to be_truthy
    end
    it 'checks if the user is an admin on the project' do
      expect(@superuser.is_admin?(:labors)).to be_truthy
      expect(@supereditor.is_admin?(:labors)).to be_falsy
      expect(@superviewer.is_admin?(:labors)).to be_falsy
      expect(@admin.is_admin?(:labors)).to be_truthy
      expect(@editor.is_admin?(:labors)).to be_falsy
      expect(@viewer.is_admin?(:labors)).to be_falsy
      expect(@guest.is_admin?(:labors)).to be_falsy
    end
    it 'checks if the user can see restricted data' do
      expect(@superuser.can_see_restricted?(:labors)).to be_falsy
      expect(@supereditor.can_see_restricted?(:labors)).to be_falsy
      expect(@superviewer.can_see_restricted?(:labors)).to be_falsy
      expect(@admin.can_see_restricted?(:labors)).to be_truthy
      expect(@editor.can_see_restricted?(:labors)).to be_truthy
      expect(@viewer.can_see_restricted?(:labors)).to be_falsy
      expect(@guest.can_see_restricted?(:labors)).to be_falsy
    end
    it 'gives a list of user projects' do
      expect(@superuser.projects).to eq(['administration'])
      expect(@supereditor.projects).to eq(['administration'])
      expect(@superviewer.projects).to eq(['administration'])
      expect(@admin.projects).to eq(['labors'])
      expect(@editor.projects).to eq(['labors'])
      expect(@viewer.projects).to eq(['labors'])
      expect(@guest.projects).to eq(['labors'])
    end
    it 'gives a list of user permissions' do
      expect(@superuser.permissions).to eq( 'administration' => { role: :admin, restricted: false } )
      expect(@supereditor.permissions).to eq( 'administration' => { role: :editor, restricted: false } )
      expect(@superviewer.permissions).to eq( 'administration' => { role: :viewer, restricted: false } )
      expect(@admin.permissions).to eq( 'labors' => { role: :admin, restricted: true } )
      expect(@editor.permissions).to eq( 'labors' => { role: :editor, restricted: true } )
      expect(@viewer.permissions).to eq( 'labors' => { role: :viewer, restricted: false } )
      expect(@guest.permissions).to eq( 'labors' => { role: :guest, restricted: false } )
    end
  end

  context "flags" do
    it 'checks flags' do
      u = Etna::User.new(
        {
          email: 'janus@two-faces.org',
          name: 'Janus Bifrons',
          perm: '',
          flags: 'doors;portals'
        },
        'xyz123randomtoken'
      )
      expect(u.has_flag?('doors')).to be_truthy
    end
    it 'checks flags if none are set' do
      u = Etna::User.new(
        {
          email: 'janus@two-faces.org',
          name: 'Janus Bifrons',
          perm: '',
        },
        'xyz123randomtoken'
      )
      expect(u.has_flag?('doors')).to be_falsy
    end
  end

  it "gives global permission to an administrator" do
    admin = Etna::User.new(
      email: 'janus@two-faces.org',
      name: 'Janus Bifrons',
      perm: 'a:administration'
    )
    expect(admin.is_admin?(:administration)).to be_truthy
    expect(admin.is_admin?(:labors)).to be_truthy
  end
end
