describe Polyphemus::TimeScanBasedEtlScanner do
  let(:scanner) do
    Polyphemus::TimeScanBasedEtlScanner.new.start_batch_state do |cursor|
      { data: data, updated_at: cursor.updated_at }
    end.execute_batch_find do |state, offset|
      state[:data].select {|r| r[:updated_at] >= Time.at(state[:updated_at].to_i) + 1}.slice(offset, limit)
    end.result_updated_at do |result|
      result[:updated_at]
    end.result_id do |result|
      result[:id]
    end
  end

  let(:limit) { 3 }
  let(:data) { [] }

  def generate_one(bump_time = false)
    @now ||= Time.now
    @now += 1 if bump_time
    @id ||= 0
    @id += 1
    data << { id: @id, updated_at: Time.at(@now.to_i + Random.rand) }
  end

  def generate_random_dataset
    data.clear
    Random.rand(1..15).times do
      generate_one(Random.rand < 0.5)
    end
  end

  def run_scan(cursor)
    all_results = []
    while true
      result = scanner.find_batch(cursor)
      all_results.push(*result)
      break if result.empty?
      yield result if block_given?
    end
    all_results
  end

  describe '#find_batch' do
    it 'never yields the same value twice' do
      40.times do
        generate_random_dataset
        results = run_scan(Polyphemus::EtlCursor.new('test-cursor'))
        expect(results).to eq(results.uniq { |r| r[:id] })
      end
    end

    it 'eventually yields all values in order' do
      40.times do
        generate_random_dataset
        results = run_scan(Polyphemus::EtlCursor.new('test-cursor'))
        expect(results).to eq(data)
      end
    end

    context 'with a persistent cursor' do
      it 'continues without duplicating, updating the cursor as it goes' do
        cursor = Polyphemus::EtlCursor.new('test-cursor')
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