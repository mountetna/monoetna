require 'date'

describe Metis::Query do

  before(:each) do
    @bucket = default_bucket('athena')

    @wisdom_folder = create_folder('athena', 'wisdom', bucket: @bucket)
    stubs.create_folder('athena', 'files', 'wisdom')

    @wisdom_file = create_file('athena', 'wisdom.txt', WISDOM, folder: @wisdom_folder)
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
    expect(results[:files].count).to eq(1)
    expect(results[:folders].count).to eq(0)

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

  it 'files and folders returned when no type specified' do
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
    expect(results[:files].count).to eq(1)
    expect(results[:folders].count).to eq(0)

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
    builder = Metis::QueryBuilder.new(
      Metis::File.where(project_name: 'athena', bucket: @bucket),
      [{
        attribute: 'updated_at',
        predicate: '>',
        value: (DateTime.now - 1).iso8601
      }])
    expect(builder.build.count).to eq(1)


    builder = Metis::QueryBuilder.new(
      Metis::File.where(project_name: 'athena', bucket: @bucket),
      [{
        attribute: 'updated_at',
        predicate: '<',
        value: (DateTime.now - 1).iso8601
      }])
    expect(builder.build.count).to eq(0)

    builder = Metis::QueryBuilder.new(
      Metis::File.where(project_name: 'athena', bucket: @bucket),
      [{
        attribute: 'updated_at',
        predicate: '<',
        value: (DateTime.now - 1).iso8601
      }, {
        attribute: 'updated_at',
        predicate: '>',
        value: (DateTime.now + 1).iso8601
      }])
    expect(builder.build.count).to eq(0)

    builder = Metis::QueryBuilder.new(
      Metis::File.where(project_name: 'athena', bucket: @bucket),
      [{
        attribute: 'updated_at',
        predicate: '<',
        value: (DateTime.now + 1).iso8601
      }, {
        attribute: 'updated_at',
        predicate: '>',
        value: (DateTime.now - 1).iso8601
      }])
    expect(builder.build.count).to eq(1)
  end

  it 'can query both string and time attributes' do
    builder = Metis::QueryBuilder.new(
      Metis::File.where(project_name: 'athena', bucket: @bucket),
      [{
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
      }])
    expect(builder.build.count).to eq(0)

    builder = Metis::QueryBuilder.new(
      Metis::File.where(project_name: 'athena', bucket: @bucket),
      [{
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
        value: 'wisd%'
      }])
    expect(builder.build.count).to eq(1)

    builder = Metis::QueryBuilder.new(
      Metis::Folder.where(project_name: 'athena', bucket: @bucket),
      [{
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
      }])
    expect(builder.build.count).to eq(0)

    builder = Metis::QueryBuilder.new(
      Metis::Folder.where(project_name: 'athena', bucket: @bucket),
      [{
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
        value: 'wisd%'
      }])
    expect(builder.build.count).to eq(1)
  end

  it 'does not return results if different type specified' do
    builder = Metis::QueryBuilder.new(
      Metis::File.where(project_name: 'athena', bucket: @bucket),
      [{
        attribute: 'updated_at',
        predicate: '<',
        value: (DateTime.now - 1).iso8601,
        type: 'folder'
      }])

    expect(builder.build.count).to eq(0)

    builder = Metis::QueryBuilder.new(
      Metis::Folder.where(project_name: 'athena', bucket: @bucket),
      [{
        attribute: 'name',
        predicate: '=~',
        value: '%publi%',
        type: 'file'
      }])
    expect(builder.build.count).to eq(0)
  end

  it 'does return results with the type flag' do
    builder = Metis::QueryBuilder.new(
      Metis::File.where(project_name: 'athena', bucket: @bucket),
      [{
        attribute: 'updated_at',
        predicate: '<',
        value: (DateTime.now + 1).iso8601,
        type: 'file'
      }])
    expect(builder.build.count).to eq(1)

    builder = Metis::QueryBuilder.new(
      Metis::Folder.where(project_name: 'athena', bucket: @bucket),
      [{
        attribute: 'name',
        predicate: '=~',
        value: '%publi%',
        type: 'folder'
      }])
    expect(builder.build.count).to eq(1)
  end
end