require 'date'

describe Metis::QueryBuilder do

  before(:each) do
    @bucket = default_bucket('athena')

    @wisdom_folder = create_folder('athena', 'wisdom', bucket: @bucket)
    stubs.create_folder('athena', 'files', 'wisdom')

    @wisdom_file = create_file('athena', 'wisdom.txt', WISDOM, folder: @wisdom_folder)
    stubs.create_file('athena', 'files', 'wisdom/wisdom.txt', WISDOM)
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

    builder = Metis::QueryBuilder.new(
      Metis::File.where(project_name: 'athena', bucket: @bucket),
      [{
        attribute: 'name',
        predicate: '=~',
        value: '/.*dOm.*/ix'
      }])
    expect(builder.build.count).to eq(1)
  end

  it 'supports querying by string equality' do
    builder = Metis::QueryBuilder.new(
      Metis::File.where(project_name: 'athena', bucket: @bucket),
      [{
        attribute: 'name',
        predicate: '=',
        value: 'dom'
      }])
    expect(builder.build.count).to eq(0)

    builder = Metis::QueryBuilder.new(
      Metis::File.where(project_name: 'athena', bucket: @bucket),
      [{
        attribute: 'name',
        predicate: '=',
        value: 'wisdom.txt'
      }])
    expect(builder.build.count).to eq(1)
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

  it 'can query by globs for folders' do
    child_folder = create_folder('athena', 'child', bucket: @bucket, folder: @wisdom_folder)
    stubs.create_folder('athena', 'files', 'wisdom/child')

    grandchild_folder = create_folder('athena', 'grandchild', bucket: @bucket, folder: child_folder)
    stubs.create_folder('athena', 'files', 'wisdom/child/grandchild')

    greatgrandchild_folder = create_folder('athena', 'greatgrandchild', bucket: @bucket, folder: grandchild_folder)
    stubs.create_folder('athena', 'files', 'wisdom/child/grandchild/greatgrandchild')

    builder = Metis::QueryBuilder.new(
      Metis::Folder.where(project_name: 'athena', bucket: @bucket),
      [{
        attribute: 'name',
        predicate: 'glob',
        value: 'child/*'
      }])
    expect(builder.build.count).to eq(2)

    builder = Metis::QueryBuilder.new(
      Metis::Folder.where(project_name: 'athena', bucket: @bucket),
      [{
        attribute: 'name',
        predicate: 'glob',
        value: 'child/gr*'
      }])
    expect(builder.build.count).to eq(1)

    builder = Metis::QueryBuilder.new(
      Metis::Folder.where(project_name: 'athena', bucket: @bucket),
      [{
        attribute: 'name',
        predicate: 'glob',
        value: 'greatgrandchild/*'
      }])
    expect(builder.build.count).to eq(0)

    builder = Metis::QueryBuilder.new(
      Metis::Folder.where(project_name: 'athena', bucket: @bucket),
      [{
        attribute: 'name',
        predicate: 'glob',
        value: 'wisdom/*'
      }])
    expect(builder.build.count).to eq(3)

    builder = Metis::QueryBuilder.new(
      Metis::Folder.where(project_name: 'athena', bucket: @bucket),
      [{
        attribute: 'name',
        predicate: 'glob',
        value: '*chil*/*'
      }])
    expect(builder.build.count).to eq(3)

    builder = Metis::QueryBuilder.new(
      Metis::Folder.where(project_name: 'athena', bucket: @bucket),
      [{
        attribute: 'name',
        predicate: 'glob',
        value: 'wisdom/**/*child*'
      }])
    expect(builder.build.count).to eq(3)
  end

  it 'can find data in folders with the same name via globbing' do
    child_folder = create_folder('athena', 'child', bucket: @bucket, folder: @wisdom_folder)
    stubs.create_folder('athena', 'files', 'wisdom/child')

    grandchild_folder = create_folder('athena', 'grandchild', bucket: @bucket, folder: child_folder)
    stubs.create_folder('athena', 'files', 'wisdom/child/grandchild')

    sibling_folder = create_folder('athena', 'sibling', bucket: @bucket, folder: @wisdom_folder)
    stubs.create_folder('athena', 'files', 'wisdom/sibling')

    grandchild2_folder = create_folder('athena', 'grandchild', bucket: @bucket, folder: sibling_folder)
    stubs.create_folder('athena', 'files', 'wisdom/sibling/grandchild')

    helmet_file = create_file('athena', 'helmet.jpg', HELMET, folder: grandchild_folder)
    stubs.create_file('athena', 'files', 'wisdom/child/grandchild/helmet.jpg', HELMET)

    shiney_helmet_file = create_file('athena', 'shiny_helmet.jpg', SHINY_HELMET, folder: grandchild2_folder)
    stubs.create_file('athena', 'files', 'wisdom/sibling/grandchild/shiny_helmet.jpg', SHINY_HELMET)

    builder = Metis::QueryBuilder.new(
      Metis::Folder.where(project_name: 'athena', bucket: @bucket),
      [{
        attribute: 'name',
        predicate: 'glob',
        value: 'grandchild'
      }])
    expect(builder.build.count).to eq(2)

    builder = Metis::QueryBuilder.new(
      Metis::Folder.where(project_name: 'athena', bucket: @bucket),
      [{
        attribute: 'name',
        predicate: 'glob',
        value: 'grandchild/*'
      }])
    expect(builder.build.count).to eq(0)

    builder = Metis::QueryBuilder.new(
      Metis::Folder.where(project_name: 'athena', bucket: @bucket),
      [{
        attribute: 'name',
        predicate: 'glob',
        value: 'chil*/*'
      }])
    expect(builder.build.count).to eq(2) # child/ + child/grandchild/

    builder = Metis::QueryBuilder.new(
      Metis::File.where(project_name: 'athena', bucket: @bucket),
      [{
        attribute: 'name',
        predicate: 'glob',
        value: 'grandchild/*'
      }])
    expect(builder.build.count).to eq(2)

    builder = Metis::QueryBuilder.new(
      Metis::File.where(project_name: 'athena', bucket: @bucket),
      [{
        attribute: 'name',
        predicate: 'glob',
        value: 'grand*/*'
      }])
    expect(builder.build.count).to eq(2)

    builder = Metis::QueryBuilder.new(
      Metis::File.where(project_name: 'athena', bucket: @bucket),
      [{
        attribute: 'name',
        predicate: 'glob',
        value: 'grandchild/*.jpg'
      }])
    expect(builder.build.count).to eq(2)

    builder = Metis::QueryBuilder.new(
      Metis::File.where(project_name: 'athena', bucket: @bucket),
      [{
        attribute: 'name',
        predicate: 'glob',
        value: 'grandchild/*.txt'
      }])
    expect(builder.build.count).to eq(0)
  end

  it 'can query by globs for files' do
    child_folder = create_folder('athena', 'child', bucket: @bucket, folder: @wisdom_folder)
    stubs.create_folder('athena', 'files', 'wisdom/child')

    grandchild_folder = create_folder('athena', 'grandchild', bucket: @bucket, folder: child_folder)
    stubs.create_folder('athena', 'files', 'wisdom/child/grandchild')

    young_wisdom_file = create_file('athena', 'young_wisdom.txt', WISDOM*2, folder: grandchild_folder)
    stubs.create_file('athena', 'files', 'wisdom/child/grandchild/young_wisdom.txt', WISDOM*2)

    builder = Metis::QueryBuilder.new(
      Metis::File.where(project_name: 'athena', bucket: @bucket),
      [{
        attribute: 'name',
        predicate: 'glob',
        value: 'wisdom/*.txt'
      }])
    expect(builder.build.count).to eq(1)
    expect(builder.build.first[:file_name]).to eq(@wisdom_file.file_name)

    builder = Metis::QueryBuilder.new(
      Metis::File.where(project_name: 'athena', bucket: @bucket),
      [{
        attribute: 'name',
        predicate: 'glob',
        value: 'wisdom/**/*.txt'
      }])
    expect(builder.build.count).to eq(2)
    expect(builder.build.first[:file_name]).to eq(@wisdom_file.file_name)
    expect(builder.build.last[:file_name]).to eq(young_wisdom_file.file_name)

    builder = Metis::QueryBuilder.new(
      Metis::File.where(project_name: 'athena', bucket: @bucket),
      [{
        attribute: 'name',
        predicate: 'glob',
        value: 'wisdom/*/*.txt'
      }])
    expect(builder.build.count).to eq(0)

    builder = Metis::QueryBuilder.new(
      Metis::File.where(project_name: 'athena', bucket: @bucket),
      [{
        attribute: 'name',
        predicate: 'glob',
        value: 'child/*.jpg'
      }])
    expect(builder.build.count).to eq(0)
  end
end