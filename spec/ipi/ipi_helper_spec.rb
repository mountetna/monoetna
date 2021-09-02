describe IpiHelper do
  before(:each) do
    copy_renaming_project
  end

  describe "rna_seq renames" do
    let(:helper) { IpiHelper.new("lib/etls/renaming/projects/test_renames.json") }

    it "returns original tube name if not in the plate renames" do
      result = helper.corrected_rna_seq_tube_name("another-random-name")
      expect(result).to eq("another-random-name")
    end

    it "returns the renamed tube name" do
      result = helper.corrected_rna_seq_tube_name("WRONG001.T1.rna.tumor")
      expect(result).to eq("RIGHT001.T1.rna.tumor")
    end
  end

  describe "rna_seq finding old 'wrong' tube_name from 'correct' tube_name" do
    let(:helper) { IpiHelper.new("lib/etls/renaming/projects/test_renames.json") }

    it "returns new tube name if not in the plate renames" do
      result = helper.incorrect_rna_seq_tube_name("another-random-name")
      expect(result).to eq("another-random-name")
    end

    it "returns the original tube name" do
      result = helper.incorrect_rna_seq_tube_name("RIGHT001.T1.rna.tumor")
      expect(result).to eq("WRONG001.T1.rna.tumor")
    end
  end
end
