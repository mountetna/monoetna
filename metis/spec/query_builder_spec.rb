require 'date'

describe Metis::QueryBuilder do

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

  it 'initializes correctly for files' do
    builder = Metis::QueryBuilder.new(
      Metis::File.where(project_name: 'athena', bucket: @bucket),
      [])
    builder.build
  end

  it 'initializes correctly for folders' do
    builder = Metis::QueryBuilder.new(
      Metis::Folder.where(project_name: 'athena', bucket: @bucket),
      [])
    builder.build
  end

  it 'supports querying by string value with sequel Like syntax' do
    builder = Metis::QueryBuilder.new(
      Metis::File.where(project_name: 'athena', bucket: @bucket),
      [{
        attribute: 'name',
        predicate: '=~',
        value: '%dom%'
      }])
    expect(builder.build.count).to eq(1)


    builder = Metis::QueryBuilder.new(
      Metis::File.where(project_name: 'athena', bucket: @bucket),
      [{
        attribute: 'name',
        predicate: '=~',
        value: '%xyz%'
      }])
    expect(builder.build.count).to eq(0)
  end

  it 'supports querying by string value with regexes' do
    builder = Metis::QueryBuilder.new(
      Metis::File.where(project_name: 'athena', bucket: @bucket),
      [{
        attribute: 'name',
        predicate: '=~',
        value: '/.*dom.*/'
      }])
    expect(builder.build.count).to eq(1)


    builder = Metis::QueryBuilder.new(
      Metis::File.where(project_name: 'athena', bucket: @bucket),
      [{
        attribute: 'name',
        predicate: '=~',
        value: '/.*xyz.*/'
      }])
    expect(builder.build.count).to eq(0)
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
end