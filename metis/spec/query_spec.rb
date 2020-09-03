require 'date'

describe Metis::Query do

  before(:each) do
    @bucket = default_bucket('athena')

    @public_folder = create_folder('athena', 'public', bucket: @bucket)
    stubs.create_folder('athena', 'files', 'public')

    @wisdom_file = create_file('athena', 'wisdom.txt', WISDOM, folder: @public_folder)
    stubs.create_file('athena', 'files', 'wisdom', 'wisdom.txt', WISDOM)
  end

  after(:each) do
    stubs.clear

    expect(stubs.contents(:athena)).to be_empty
  end

  it 'initializes correctly' do
    query = Metis::Query.new(
      project_name: 'athena',
      bucket: @bucket,
      params: []
    )
  end

  it 'supports querying files by string value' do
    query = Metis::Query.new(
      project_name: 'athena',
      bucket: @bucket,
      params: [{
        attribute: 'name',
        predicate: '=~',
        value: '%dom%',
        type: 'file'
      }]
    )

    results = query.execute
    expect(results[:files].count).to eq(1)
    expect(results[:folders].count).to eq(0)

    query = Metis::Query.new(
      project_name: 'athena',
      bucket: @bucket,
      params: [{
        attribute: 'name',
        predicate: '=~',
        value: '%xyz%',
        type: 'file'
      }]
    )
    results = query.execute
    expect(results[:files].count).to eq(0)
    expect(results[:folders].count).to eq(0)
  end

  it 'supports querying folders by string value' do
    query = Metis::Query.new(
      project_name: 'athena',
      bucket: @bucket,
      params: [{
        attribute: 'name',
        predicate: '=~',
        value: '%pub%',
        type: 'folder'
      }]
    )

    results = query.execute
    expect(results[:files].count).to eq(0)
    expect(results[:folders].count).to eq(1)

    query = Metis::Query.new(
      project_name: 'athena',
      bucket: @bucket,
      params: [{
        attribute: 'name',
        predicate: '=~',
        value: '%priv%',
        type: 'folder'
      }]
    )

    results = query.execute
    expect(results[:files].count).to eq(0)
    expect(results[:folders].count).to eq(0)
  end

  it 'both files and folders returned when no type specified' do
    query = Metis::Query.new(
      project_name: 'athena',
      bucket: @bucket,
      params: [{
        attribute: 'name',
        predicate: '=~',
        value: '%pub%'
      }]
    )

    results = query.execute
    expect(results[:files].count).to eq(0)
    expect(results[:folders].count).to eq(1)

    query = Metis::Query.new(
      project_name: 'athena',
      bucket: @bucket,
      params: [{
        attribute: 'name',
        predicate: '=~',
        value: '%i%'
      }]
    )

    results = query.execute
    expect(results[:files].count).to eq(1)
    expect(results[:folders].count).to eq(1)

    query = Metis::Query.new(
      project_name: 'athena',
      bucket: @bucket,
      params: [{
        attribute: 'name',
        predicate: '=~',
        value: '%w%'
      }]
    )

    results = query.execute
    expect(results[:files].count).to eq(1)
    expect(results[:folders].count).to eq(0)
  end

  it 'supports querying by time value and range' do
    query = Metis::Query.new(
      project_name: 'athena',
      bucket: @bucket,
      params: [{
        attribute: 'updated_at',
        predicate: '>',
        value: (DateTime.now - 1).iso8601
      }]
    )
    results = query.execute
    expect(results[:files].count).to eq(1)
    expect(results[:folders].count).to eq(1)

    query = Metis::Query.new(
      project_name: 'athena',
      bucket: @bucket,
      params: [{
        attribute: 'updated_at',
        predicate: '<',
        value: (DateTime.now - 1).iso8601
      }]
    )
    results = query.execute
    expect(results[:files].count).to eq(0)
    expect(results[:folders].count).to eq(0)

    query = Metis::Query.new(
      project_name: 'athena',
      bucket: @bucket,
      params: [{
        attribute: 'updated_at',
        predicate: '<',
        value: (DateTime.now - 1).iso8601
      }, {
        attribute: 'updated_at',
        predicate: '>',
        value: (DateTime.now + 1).iso8601
      }]
    )
    results = query.execute
    expect(results[:files].count).to eq(0)
    expect(results[:folders].count).to eq(0)

    query = Metis::Query.new(
      project_name: 'athena',
      bucket: @bucket,
      params: [{
        attribute: 'updated_at',
        predicate: '>',
        value: (DateTime.now - 1).iso8601
      }, {
        attribute: 'updated_at',
        predicate: '<',
        value: (DateTime.now + 1).iso8601
      }]
    )
    results = query.execute
    expect(results[:files].count).to eq(1)
    expect(results[:folders].count).to eq(1)
  end

  it 'can query both string and time attributes' do
    query = Metis::Query.new(
      project_name: 'athena',
      bucket: @bucket,
      params: [{
        attribute: 'updated_at',
        predicate: '<',
        value: (DateTime.now + 1).iso8601
      }, {
        attribute: 'updated_at',
        predicate: '>',
        value: (DateTime.now - 1).iso8601
      }, {
        attribute: 'name',
        predicate: '=~',
        value: '%foo%'
      }]
    )
    results = query.execute
    expect(results[:files].count).to eq(0)
    expect(results[:folders].count).to eq(0)

    query = Metis::Query.new(
      project_name: 'athena',
      bucket: @bucket,
      params: [{
        attribute: 'updated_at',
        predicate: '<',
        value: (DateTime.now + 1).iso8601
      }, {
        attribute: 'updated_at',
        predicate: '>',
        value: (DateTime.now - 1).iso8601
      }, {
        attribute: 'name',
        predicate: '=~',
        value: '%wis%'
      }]
    )
    results = query.execute
    expect(results[:files].count).to eq(1)
    expect(results[:folders].count).to eq(0)

    query = Metis::Query.new(
      project_name: 'athena',
      bucket: @bucket,
      params: [{
        attribute: 'updated_at',
        predicate: '<',
        value: (DateTime.now + 1).iso8601
      }, {
        attribute: 'updated_at',
        predicate: '>',
        value: (DateTime.now - 1).iso8601
      }, {
        attribute: 'name',
        predicate: '=~',
        value: '%wis%',
        type: 'folder'
      }]
    )
    results = query.execute
    expect(results[:files].count).to eq(1)
    expect(results[:folders].count).to eq(0)
  end
end