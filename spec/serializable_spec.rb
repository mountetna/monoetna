require_relative '../lib/etna/multipart_serializable_nested_hash'
require_relative '../lib/etna/json_serializable_struct'

describe Etna::MultipartSerializableNestedHash do
  class MyHash < Hash
    include Etna::MultipartSerializableNestedHash
  end

  def generate_value
    [
        -> { generate_scalar },
        -> { generate_array },
        -> { generate_hash },
    ].sample.call
  end

  def generate_scalar
    [
        "100",
        "movie",
        1.0,
        "blah blah blah ho&ho&di",
        StringIO.new("test-string=io"),
    ].sample
  end

  def generate_array(depth = 0)
    if depth > 4
      return generate_scalar
    end

    inner = [
        -> { generate_scalar },
        -> { generate_hash(depth + 1) },
    ].sample
    ([nil] * Random.rand(1..5)).map { inner.call }
  end

  def generate_hash(depth = 0)
    if depth > 4
      return generate_scalar
    end

    inner = [
        -> { generate_scalar },
        -> { generate_array(depth + 1) },
        -> { generate_hash(depth + 1) },
    ].sample
    ([nil] * Random.rand(1..5)).map { [generate_key, inner.call] }.to_h
  end

  def generate_key
    [
        "a".."z",
        "aa".."zz",
        "aaa"..."zzz"
    ].sample.to_a.sample
  end

  def get_expected_string_version(v)
    if v.is_a? Hash
      v.map { |k, v| [k, get_expected_string_version(v)] }.to_h
    elsif v.is_a? Array
      v.each_with_index.map { |v, i| get_expected_string_version(v).tap { |v| v["_idx"] = i.to_s if v.is_a?(Hash) } }
    elsif v.respond_to? :read
      v.rewind
      v.read
    else
      v.to_s
    end
  end

  def get_unpacked_result(v)
    if v.is_a? Hash
      if v.include? :tempfile
        v[:tempfile].read
      else
        v.map { |k, v| [k, get_unpacked_result(v)] }.to_h
      end
    elsif v.is_a? Array
      v.map { |v| get_unpacked_result(v) }
    else
      v
    end
  end

  describe 'conversion of values' do
    it 'should be symmetrical with multipart parsing' do
      100.times do
        base_value = generate_value

        hash = MyHash.new
        hash.update(inner: base_value)

        content = hash.encode_multipart_content
        multipart = Net::HTTP::Post::Multipart.new "https://test", content
        io = multipart.body_stream
        content_length = io.read.length
        io.rewind


        tempfile = Rack::Multipart::Parser::TEMPFILE_FACTORY
        bufsize = Rack::Multipart::Parser::BUFSIZE
        params = Rack::Utils.default_query_parser
        info = Rack::Multipart::Parser.parse io, content_length, multipart['content-type'], tempfile, bufsize, params

        expect(get_unpacked_result(info.params)).to eq({"inner" => get_expected_string_version(base_value)})
      end
    end
  end
end

describe Etna::JsonSerializableStruct do
  class MyStruct < Struct.new(:a, :b, keyword_init: true)
    include Etna::JsonSerializableStruct
  end

  it 'deeply serializes' do
    expect(JSON.parse(MyStruct.new(
        a: [1, 2, {b: 3}, [MyStruct.new(b: {c: 4})]],
        b: "value"
    ).to_json)).to eq({"a" => [1, 2, {"b" => 3}, [{"b" => {"c" => 4}}]], "b" => "value"})
  end
end