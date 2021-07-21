require "securerandom"

describe Polyphemus::FilenameScanBasedEtlScanner do
  class TestRsyncChange
    def initialize(filename)
      @filename = filename
    end

    def filename
      @filename
    end
  end

  let(:simple_scanner) do
    Polyphemus::FilenameScanBasedEtlScanner.new.execute_batch_find do
      data
    end.result_id do |result|
      result.filename
    end.result_updated_at do |result|
      "new"
    end
  end

  let(:limit) { 3 }
  let(:data) { [] }

  def generate_one
    data << TestRsyncChange.new(SecureRandom.hex)
  end

  def generate_random_dataset
    data.clear
    Random.rand(0..15).times do
      generate_one
    end
  end

  def run_scan(cursor)
    all_results = []
    while true
      result = simple_scanner.find_batch(cursor)
      all_results.push(*result)
      break if result.empty?
      yield result if block_given?
    end
    all_results
  end

  describe "#find_batch" do
    it "never yields the same value twice" do
      40.times do
        generate_random_dataset
        results = run_scan(Polyphemus::EtlCursor.new("test-cursor"))
        expect(results).to eq(results.uniq { |r| r.filename })
      end
    end

    it "eventually yields all values in order" do
      40.times do
        generate_random_dataset
        results = run_scan(Polyphemus::EtlCursor.new("test-cursor"))
        expect(results).to eq(data)
      end
    end

    it "finds updated records only if file name changed" do
      data.clear

      updated_entry = TestRsyncChange.new(SecureRandom.hex)

      data << updated_entry

      static_hash = SecureRandom.hex
      static_entry = TestRsyncChange.new(static_hash)

      data << static_entry

      cursor = Polyphemus::EtlCursor.new(
        "test-cursor",
        { "seen_ids" => [[SecureRandom.hex, "new"], [static_hash, "new"]] }
      )

      results = run_scan(cursor)
      expect(results).to match_array([updated_entry, static_entry])

      # Non infinite, the cursor is updated
      results = run_scan(cursor)
      expect(results).to eq([])
    end

    context "with a persistent cursor" do
      it "continues without duplicating, updating the cursor as it goes" do
        cursor = Polyphemus::EtlCursor.new("test-cursor")
        all_time_results = []

        40.times do
          generate_random_dataset
          results = run_scan(cursor)
          all_time_results.push(*results)
          expect(results).to eq(data)
        end

        expect(all_time_results).to eq(all_time_results.uniq { |r| r.filename })
      end
    end
  end
end
