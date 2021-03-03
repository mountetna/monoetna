describe Etna::Filesystem do
  describe Etna::Filesystem::GeneAsperaCliFilesystem do
    let(:filesystem) do
      Etna::Filesystem::GeneAsperaCliFilesystem.new(
          ascli_bin: `which ascli`.chomp,
          ascp_bin: "/home/aspera/.aspera/connect/bin/ascp",
          host: "demo.asperasoft.com",
          username: "aspera",
          password: "demoaspera",
      )
    end

    xit "works, sort of" do
      inner_dir = "/Upload/some/inner/dir"
      size = 1024 * 1024

      filesystem.mkdir_p(inner_dir)

      filesystem.with_writeable(File.join(inner_dir, "test.txt"), size_hint: size) do |file|
        4.times do |i|
          file.write("z" * (size / 4))
          sleep 1
        end
      end

      filesystem.with_readable(File.join(inner_dir, "test.txt")) do |file|
        expect(file.read.length).to eql(size)
      end
    end
  end

  [
      Etna::Filesystem.new,
      Etna::Filesystem::Mock.new do |fname, opts|
        StringIO.new
      end,
      Etna::Filesystem::AsperaCliFilesystem.new(
          ascli_bin: `which ascli`.chomp,
          ascp_bin: "/home/aspera/.aspera/connect/bin/ascp",
          host: "demo.asperasoft.com",
          username: "aspera",
          password: "demoaspera",
      )
  ].each do |fs|
    describe "#{fs.class.name}" do
      xit "works" do
        dir = fs.tmpdir
        begin
          inner_dir = File.join(dir, "some/inner/dir")
          expect(fs.exist?(inner_dir)).to eq(false)
          fs.mkdir_p(inner_dir)
          expect(fs.exist?(dir)).to eq(true)
          expect(fs.exist?(inner_dir)).to eq(true)
          expect(fs.exist?("/never/a/directory")).to eq(false)

          size = 1024 * 1024

          fs.with_writeable(File.join(inner_dir, "test.txt"), size_hint: size) do |file|
            4.times do |i|
              file.write("z" * (size / 4))
              sleep 1
            end
          end

          fs.with_readable(File.join(inner_dir, "test.txt")) do |file|
            expect(file.read.length).to eql(size)
          end

          expect do
            fs.with_writeable(File.join(inner_dir, "vvv.txt"), size_hint: 100) do |file|
              raise "Oh no!"
            end
          end.to raise_error

          expect do
            fs.with_readable(File.join(inner_dir, "test.txt")) do |file|
              raise "Oh no!"
            end
          end.to raise_error
        ensure
          fs.rm_rf(dir)
        end
      end
    end
  end
end