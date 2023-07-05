describe Magma::Loader do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  before(:each) do
    @user = Etna::User.new(
        email: 'zeus@mountolympus.org',
        name: 'Zeus',
        perm: 'A:labors'
    )
  end

  it 'bulk-creates records' do
    loader = Magma::Loader.new(@user, 'labors')
    loader.push_record(Labors::Labor, 'Nemean Lion', number: 1, completed: true)
    loader.push_record(Labors::Labor, 'Lernean Hydra', number: 2, completed: false)
    loader.push_record(Labors::Labor, 'Augean Stables', number: 5, completed: false)

    loader.dispatch_record_set

    expect(Labors::Labor.count).to eq(3)
  end

  it 'does not bulk-creates records when dry_run flag set' do
    loader = Magma::Loader.new(@user, 'labors', dry_run: true)
    loader.push_record(Labors::Labor, 'Nemean Lion', number: 1, completed: true)
    loader.push_record(Labors::Labor, 'Lernean Hydra', number: 2, completed: false)
    loader.push_record(Labors::Labor, 'Augean Stables', number: 5, completed: false)

    payload = loader.dispatch_record_set

    expect(Labors::Labor.count).to eq(0)
    expect(payload[:models][:labor][:documents].keys.length).to eq(3)
  end

  it 'bulk-creates records with different data' do
    loader = Magma::Loader.new(@user, 'labors')
    loader.push_record(Labors::Labor, 'Nemean Lion', number: 1)
    loader.push_record(Labors::Labor, 'Lernean Hydra', completed: false)
    loader.push_record(Labors::Labor, 'Augean Stables', number: 5)

    loader.dispatch_record_set

    expect(Labors::Labor.count).to eq(3)

    lion = Labors::Labor[name: 'Nemean Lion']
    hydra = Labors::Labor[name: 'Lernean Hydra']

    expect(lion.number).to eq(1)
    expect(lion.completed).to be_nil
    expect(hydra.number).to be_nil
    expect(hydra.completed).to eq(false)
  end

  it 'bulk-updates records' do
    lion = create(:labor, name: 'Nemean Lion', number: 1, completed: false)
    hydra = create(:labor, name: "Lernean Hydra", number: 2, completed: false)

    loader = Magma::Loader.new(@user, 'labors')
    loader.push_record(Labors::Labor, 'Nemean Lion', number: 1, completed: true)
    loader.push_record(Labors::Labor, 'Augean Stables', number: 5, completed: false)

    loader.dispatch_record_set
    lion.refresh

    expect(Labors::Labor.count).to eq(3)
    expect(lion.completed).to eq(true)
  end

  it 'does not bulk-updates records when dry_run flag set' do
    lion = create(:labor, name: 'Nemean Lion', number: 1, completed: false)
    hydra = create(:labor, name: "Lernean Hydra", number: 2, completed: false)

    expect(Labors::Labor.count).to eq(2)

    loader = Magma::Loader.new(@user, 'labors', dry_run: true)
    loader.push_record(Labors::Labor, 'Nemean Lion', number: 1, completed: true)
    loader.push_record(Labors::Labor, 'Augean Stables', number: 5, completed: false)

    payload = loader.dispatch_record_set
    lion.refresh

    expect(Labors::Labor.count).to eq(2)
    expect(lion.completed).to eq(false)
    expect(payload[:models][:labor][:documents][lion.name][:completed]).to eq(true)
  end

  it 'bulk-updates records with different data' do
    lion = create(:labor, name: 'Nemean Lion', number: 1, completed: false)
    hydra = create(:labor, name: "Lernean Hydra", number: 3, completed: false)

    loader = Magma::Loader.new(@user, 'labors')
    loader.push_record(Labors::Labor, 'Nemean Lion', completed: true)
    loader.push_record(Labors::Labor, 'Lernean Hydra', number: 2)

    loader.dispatch_record_set

    expect(Labors::Labor.count).to eq(2)

    lion = Labors::Labor[name: 'Nemean Lion']
    hydra = Labors::Labor[name: 'Lernean Hydra']

    expect(lion.number).to eq(1)
    expect(lion.completed).to eq(true)
    expect(hydra.number).to eq(2)
    expect(hydra.completed).to eq(false)
  end

  it 'validates records' do
    loader = Magma::Loader.new(@user, 'labors')
    loader.push_record(Labors::Monster, 'Nemean Lion', species: 'Lion')

    expect { loader.dispatch_record_set }.to raise_error(Magma::LoadFailed)
  end

  it 'creates associations' do
    loader = Magma::Loader.new(@user, 'labors')
    loader.push_record(Labors::Labor, 'Nemean Lion', species: 'Lion', prize: [ '::temphide' ]) 
    loader.push_record(Labors::Prize, '::temphide', labor: 'Nemean Lion', name: 'hide')

    loader.dispatch_record_set
    lion = Labors::Labor[name: 'Nemean Lion']

    expect(lion.prize.first.name).to eq('hide')
  end

  it 'date-shifts payload only, during a dry-run' do
    set_date_shift_root('monster', true)
    lion = create(:labor, name: 'Nemean Lion', number: 1, completed: false)
    lion_monster = create(:monster, name: 'Nemean Lion', labor: lion)
    
    john_doe = create(:victim, name: "John Doe", monster: lion_monster)

    expect(john_doe[:birthday]).to eq(nil)
    expect(Labors::Victim.count).to eq(1)

    loader = Magma::Loader.new(@user, 'labors', dry_run: true)
    loader.push_record(Labors::Victim, 'John Doe', birthday: '2000-01-01')
    loader.push_record(Labors::Victim, 'Jane Doe', monster: lion_monster.name, birthday: '2000-01-01')

    payload = loader.dispatch_record_set
    john_doe.refresh

    expect(Labors::Victim.count).to eq(1)
    expect(john_doe[:birthday]).to eq(nil)
    expect(payload[:models][:victim][:documents][john_doe.name][:birthday].iso8601).not_to eq(iso_date_str('2000-01-01'))
    expect(payload[:models][:victim][:documents]['Jane Doe'][:birthday].iso8601).not_to eq(iso_date_str('2000-01-01'))

    set_date_shift_root('monster', false)
  end

  it 'can find parent models when hierarchical grammar exists' do
    grammar = create(:grammar, { project_name: 'labors', version_number: 1, config: HIERARCHY_GRAMMAR_CONFIG, comment: 'update' })
    project_identifier = create_identifier("The Twelve Labors of Hercules", rule: 'project', grammar: grammar)
    labor_identifier = create_identifier("The Nemean Lion", rule: 'labor', grammar: grammar)
    monster_identifier = create_identifier("Neamon Lion", rule: 'monster', grammar: grammar)
    victim_identifier = create_identifier("LABORS-LION-NEAMON-H2-C1", rule: 'victim', grammar: grammar)

    loader = Magma::Loader.new(@user, 'labors')
    parent_models = loader.find_parent_models(victim_identifier.identifier)

    expect(parent_models.length).to eq(4)

    expect(parent_models[0][:model]).to eq(:project)
    expect(parent_models[0][:identifier]).to eq(project_identifier.identifier)
    expect(parent_models[0][:parent_model]).to eq(nil)

    expect(parent_models[1][:model]).to eq(:labor)
    expect(parent_models[1][:identifier]).to eq(labor_identifier.identifier)
    expect(parent_models[1][:parent_model]).to eq(:project)

    binding.pry
    expect(parent_models[2][:model]).to eq(:monster)
    expect(parent_models[2][:identifier]).to eq(monster_identifier.identifier)
    expect(parent_models[2][:parent_model]).to eq(:labor)

    expect(parent_models[3][:model]).to eq(:victim)
    expect(parent_models[3][:identifier]).to eq(victim_identifier.identifier)
    expect(parent_models[3][:parent_model]).to eq(:monster)

  end

  it 'cannot find parent models when the grammar is non hierarchical' do
    # The grammar does not describe a MONSTER token
    grammar = create(:grammar, { project_name: 'labors', version_number: 1, config: VALID_GRAMMAR_CONFIG, comment: 'update' })

    create_identifier("Nemean Lion", rule: 'monster', grammar: grammar)
    victim_identifier = create_identifier("LABORS-LION-H2-C1", rule: 'victim', grammar: grammar)

    loader = Magma::Loader.new(@user, 'labors')
    parent_models = loader.find_parent_models(victim_identifier.identifier)
    expect(parent_models.empty?).to eq(true)
  end


end
