require 'sequel'

describe Polyphemus::EtlCursor do
  describe '#load_from_db' do
    let(:name) { 'cursor_a' }
    let(:cursor) { Polyphemus::EtlCursor.new(name) }

    context 'when no existing value by its name exists' do
      it 'retains any local data' do
        time = Time.now
        cursor.updated_at = time
        cursor['abc'] = 1
        expect(cursor.value).to eql({ 'abc' => 1 })
        expect(cursor.updated_at).to eql(time)

        cursor.load_from_db

        expect(cursor.value).to eql({ 'abc' => 1})
        expect(cursor.updated_at).to eql(time)
      end
    end

    context 'when an existing value by the cursor name exists' do
      before(:each) do
        existing = Polyphemus::EtlCursor.new(name)
        @time = Time.now - 1000
        existing.updated_at = @time
        existing[:some_key] = 123

        existing.save_to_db
        cursor.load_from_db
      end

      it 'replaces with values from the existing cursor' do
        expect(cursor.value).to eql({'some_key' => 123})
        # Unfortunately, the psql serialization reduces fidelity here by dropping sec_fractions.
        expect(cursor.updated_at).to be_between(@time - 1, @time)
      end
    end
  end

  describe 'saving to db' do
    before(:each) do
      @cursor_a = Polyphemus::EtlCursor.new('name_a')
      @cursor_b = Polyphemus::EtlCursor.new('name_b')

      @cursor_a.updated_at = Time.now
      @cursor_b.updated_at = Time.now

      @cursor_a[:v] = 1
      @cursor_b[:c] = 2

      @cursor_a.save_to_db
      @cursor_b.save_to_db
    end

    it 'inserts by name' do
      expect(Polyphemus::EtlCursor.new('name_a').load_from_db.value['v']).to eql(1)
      expect(Polyphemus::EtlCursor.new('name_b').load_from_db.value['c']).to eql(2)
    end

    context 'over an existing entry' do
      before(:each) do
        @cursor_a[:d] = 4
        @time = Time.now - 500
        @cursor_a.updated_at = @time
        @cursor_a.save_to_db
      end

      it 'updates by name' do
        expect(Polyphemus::EtlCursor.new('name_a').load_from_db.value['d']).to eql(4)
        expect(Polyphemus::EtlCursor.new('name_a').load_from_db.updated_at).to be_between(@time - 1, @time)
        expect(Polyphemus::EtlCursor.new('name_b').load_from_db.value['c']).to eql(2)
      end

      it 'uses optimistic locking' do
        expect(Polyphemus::EtlCursor.new('name_a').load_from_db.value['d']).to eql(4)
        c = Polyphemus::EtlCursor.new('name_a').load_from_db
        c[:d] = 10
        c.save_to_db
        @cursor_a[:d] = 9
        expect(@cursor_a.save_to_db).to be_falsey
        expect(@cursor_a[:d]).to eq(10)
        expect(Polyphemus::EtlCursor.new('name_a').load_from_db.value['d']).to eql(10)
      end
    end
  end
end

describe Polyphemus::EtlCursorGroup do
  let(:group) { Polyphemus::EtlCursorGroup.new(cursors) }
  let(:cursors) { [] }
  describe '#with_next' do
    it 'does not yield when no cursors are given' do
      @yielded = false
      group.with_next { @yielded = true }
      expect(@yielded).to be_falsey
    end

    context 'with some cursors' do
      let(:cursors) do
        [
            Polyphemus::EtlCursor.new('a').tap { |c| c.updated_at = Time.now - 5 },
            Polyphemus::EtlCursor.new('b').tap { |c| c.updated_at = Time.now - 3 },
            Polyphemus::EtlCursor.new('c').tap { |c| c.updated_at = Time.now - 9 },
        ]
      end

      it 'yields the oldest cursor first' do
        @yielded = false
        expect(group.with_next { |v| @yielded = v; true }).to eql(true)
        expect(@yielded).to eql(cursors.last)
      end

      it 'continues yielding on false' do
        @yielded = []
        expect(group.with_next { |v| @yielded << v; false }).to eql(false)
        expect(@yielded).to eql([cursors[2], cursors[0], cursors[1]])
      end
    end
  end
end
