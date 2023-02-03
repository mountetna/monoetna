describe Etna::Route do
  include Rack::Test::Methods

  it "correctly escapes filename for a route" do
    file_route = Etna::Route.new(
      'GET',
      '/files/:file_name',
      {}
    )

    path = file_route.path({
      file_name: 'test file+symbols.txt'
    })
    expect(path).to eq(
      "/files/test%20file%2Bsymbols.txt"
    )
  end
end
