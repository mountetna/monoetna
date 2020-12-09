describe Etna::CsvImporter do
  describe '#each_csv_row' do
    it 'maps both string headings and integer indexable headings to each row' do
      io = StringIO.new(<<-STRING)
,heading1,heading2,heading3
1,2,3,4,5,6,7,8
      STRING

      Etna::CsvImporter.new.each_csv_row(input_io: io) do |row|
        expect(row.keys.map(&:inspect)).to eql(["0", ":heading1", ":heading2", ":heading3", "4", "5", "6", "7"])
      end
    end
  end

  describe '#feed_rows_into' do
    def generator(result)
      p = yield
      result << (p.row[4] + "#{p.row[:lineno]}")
      p = yield
      result << (p.row[5] + "#{p.row[:lineno]}")
    end

    it 'allows for generator based processing backed by each_csv_row' do
      io = StringIO.new(<<-STRING)
,heading1,heading2,heading3
1,2,3,4,5,6,7,8
      STRING

      expect([].tap do |result|
        Etna::CsvImporter.new.feed_rows_into(enum_for(:generator, result), input_io: io) do |row, lineno|
          row[:lineno] = lineno
        end
      end).to eql(["52"])
    end
  end
end