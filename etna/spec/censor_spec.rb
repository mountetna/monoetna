describe Etna::Censor do
  it "redacts keys if they are hashes" do
    censor = Etna::Censor.new([:secret_hash])

    result = censor.redact(:secret_hash, {
      something: ["I don't want to reveal"],
    })

    expect(result).to eq("*")
  end

  it "redacts keys if they are strings" do
    censor = Etna::Censor.new([:secret_value])

    result = censor.redact(:secret_value, "I don't want to reveal")

    expect(result).to eq("*")
  end
end
