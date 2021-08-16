describe IpiHelper do
  describe "rna_seq renames" do
    let(:helper) { IpiHelper.new }

    it "returns original tube name if plate not included" do
      result = helper.corrected_rna_seq_tube_name("PlateX", "some-random-tube-name")
      expect(result).to eq("some-random-tube-name")
    end

    it "returns original tube name if not in the plate renames" do
      result = helper.corrected_rna_seq_tube_name("Plate4", "another-random-name")
      expect(result).to eq("another-random-name")
    end

    it "returns the renamed tube name" do
      result = helper.corrected_rna_seq_tube_name("Plate10", "IPIMEL271.T1.rna.tumor")
      expect(result).to eq("IPIMEL273.T1.rna.tumor")
    end
  end

  describe "rna_seq finding old 'wrong' tube_name from 'correct' tube_name" do
    let(:helper) { IpiHelper.new }

    it "returns original tube name if plate not included" do
      result = helper.incorrect_rna_seq_tube_name("PlateX", "some-random-tube-name")
      expect(result).to eq("some-random-tube-name")
    end

    it "returns original tube name if not in the plate renames" do
      result = helper.incorrect_rna_seq_tube_name("Plate4", "another-random-name")
      expect(result).to eq("another-random-name")
    end

    it "returns the renamed tube name" do
      result = helper.incorrect_rna_seq_tube_name("Plate10", "IPIMEL273.T1.rna.tumor")
      expect(result).to eq("IPIMEL271.T1.rna.tumor")
    end
  end
end
