require "securerandom"

describe Polyphemus::HashScanBasedEtlScanner do
  let(:random_scanner) do
    Polyphemus::HashScanBasedEtlScanner.new.start_batch_state do |cursor|
      { data: data, updated_at: cursor.updated_at }
    end.execute_batch_find do |state, expansion_number|
      data = state[:data]
      # Simulate changes to the underlying data
      if concurrent_updates
        Random.rand(0..3).times do
          next if data.empty?
          idx = Random.rand(0...(data.length))
          update = data[idx].dup
          update[:updated_at] = (data.map { |r| r[:updated_at] }.max + 1)
          update[:hash] = SecureRandom.hex
          concurrent_updates << update
        end
      end
      data.sort { |a, b| a[:updated_at] <=> b[:updated_at] }
      data.select { |r| r[:updated_at] >= Time.at(state[:updated_at].to_i) }.slice(0, limit * expansion_number)
    end.result_updated_at do |result|
      result[:updated_at]
    end.result_id do |result|
      result[:id]
    end.result_file_hashes do |result|
      result[:hash]
    end
  end

  let(:simple_scanner) do
    Polyphemus::HashScanBasedEtlScanner.new.start_batch_state do |cursor|
      { data: data, updated_at: cursor.updated_at }
    end.execute_batch_find do |state, expansion_number|
      data = state[:data]
      data.sort { |a, b| a[:updated_at] <=> b[:updated_at] }
      data.select { |r| r[:updated_at] >= Time.at(state[:updated_at].to_i) }.slice(0, limit * expansion_number)
    end.result_updated_at do |result|
      result[:updated_at]
    end.result_id do |result|
      result[:id]
    end.result_file_hashes do |result|
      result[:hash]
    end
  end

  let(:concurrent_updates) { nil }
  let(:limit) { 3 }
  let(:data) { [] }

  def generate_one(bump_time = false)
    @now ||= Time.now
    @now += 1 if bump_time
    @id ||= 0
    @id += 1
    data << { id: @id, updated_at: Time.at(@now.to_i + Random.rand), hash: SecureRandom.hex }
  end

  def generate_random_dataset
    data.clear
    Random.rand(0..15).times do
      generate_one(Random.rand < 0.5)
    end
  end

  def run_scan(cursor, random = true)
    all_results = []
    while true
      result = (random ? random_scanner : simple_scanner).find_batch(cursor)
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
        expect(results).to eq(results.uniq { |r| r[:id] })
      end
    end

    it "eventually yields all values in order" do
      40.times do
        generate_random_dataset
        results = run_scan(Polyphemus::EtlCursor.new("test-cursor"))
        expect(results).to eq(data)
      end
    end

    it "finds updated records only if hash changed" do
      now = Time.now
      Timecop.freeze
      data.clear

      updated_entry = {
        id: 1,
        updated_at: now - 10,
        hash: SecureRandom.hex,
      }

      data << updated_entry

      static_hash = SecureRandom.hex
      data << {
        id: 2,
        updated_at: now,
        hash: static_hash,
      }

      cursor = Polyphemus::EtlCursor.new(
        "test-cursor",
        Time.now - 100,
        { "seen_ids" => [[1, SecureRandom.hex], [2, static_hash]] }
      )

      results = run_scan(cursor, false)
      expect(results).to eq([updated_entry])

      # Non infinite, the cursor is updated
      results = run_scan(cursor, false)
      expect(results).to eq([])

      Timecop.return
    end

    context "with concurrent updates" do
      let(:concurrent_updates) { [] }

      it "eventually yields all values" do
        40.times do
          generate_random_dataset
          results = run_scan(Polyphemus::EtlCursor.new("test-cursor"))

          result_ids = results.map { |r| r[:id] }
          data.each { |r| result_ids.include?(r[:id]) }
        end
      end
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

        expect(all_time_results).to eq(all_time_results.uniq { |r| r[:id] })
      end
    end
  end
end
