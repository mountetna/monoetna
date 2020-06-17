require_relative '../lib/commands'

# A note on Etna::Command#setup
#
# Unfortunately Rack already sets up the application - this may cause us
# headache some day, since the commands are meant to call setup(config)
#
# If it ever does, we should figure out how to exercise setup here.
#
# For now we skip setup

describe Janus::Console do
  it 'starts a console' do
    require 'irb'
    allow(IRB).to receive(:start)

    c = Janus::Console.new

    c.execute

    expect(IRB).to have_received(:start)
  end
end

describe Janus::AddProject do
  it 'adds a project to the database' do
    project_name = 'augean_stables'
    project_name_full = 'Augean Stables'

    command = Janus::AddProject.new
    command.execute(project_name, project_name_full)

    project = Project.first

    expect(project.project_name).to eq(project_name)
  end

  it 'updates if the project exists' do
    project_name = 'augean_stables'
    project_name_full = 'Augean Stables'
    project = create(:project, project_name: project_name, project_name_full: 'Augaian Stables')

    command = Janus::AddProject.new
    command.execute(project_name, project_name_full)

    project.refresh

    expect(project.project_name_full).to eq(project_name_full)
  end
end
describe Janus::AddUser do
  it 'adds a user to the database' do
    email = 'test_user@test.edu'

    command = Janus::AddUser.new
    command.execute(email, 'Janus', 'Two-faces')

    user = User.first

    expect(user.email).to eq(email)
  end

  it 'updates if the user exists' do
    email = 'test_user@test.edu'
    last_name = 'Two-faces'
    user = create(:user, email: email, first_name: 'Janus', last_name: 'One-face')

    command = Janus::AddUser.new
    command.execute(email, 'Janus', last_name)

    user.refresh

    expect(user.last_name).to eq(last_name)
  end
end

describe Janus::Permit do
  it 'gives a user a specific permission on a project' do
    email = 'hercules@test.edu'
    project_name = 'augean_stables'
    role = 'viewer'

    user = create(:user, email: email)
    project = create(:project, project_name: project_name, project_name_full: 'Augean Stables')

    command = Janus::Permit.new
    command.execute(email, project_name, role)

    perm = Permission.first

    expect(perm.user.email).to eq(email)
    expect(perm.project.project_name).to eq(project_name)
    expect(perm.role).to eq(role)
  end
end
