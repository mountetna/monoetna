require_relative '../lib/magma'
require 'yaml'

describe Magma::Attribute do
  describe "#json_template" do
    it "includes attribute defaults" do
      attribute = Magma::Attribute.new(
        attribute_name: "name",
        format_hint: "Hint"
      )

      template = attribute.json_template

      expect(template[:display_name]).to eq("Name")
      expect(template[:format_hint]).to eq("Hint")
    end

    it "includes validation arrays as options" do
      attribute = Magma::Attribute.new(
        attribute_name: "fruits",
        validation: { type: "Array", value: ["apple", "banana"] }
      )

      template = attribute.json_template

      expect(template[:options]).to eq(["apple", "banana"])
    end

    it "includes validation regexes as match" do
      attribute = Magma::Attribute.new(
        attribute_name: "fruits",
        validation: { type: "Regexp", value: /.*/ }
      )

      template = attribute.json_template

      expect(template[:match]).to eq(".*")
    end

    it "contains a ValidationOjbect" do
      attribute = Magma::Attribute.new(
        attribute_name: "name",
        validation: { type: "Regexp", value: /^[a-zA-Z]{1}$/ }
      )

      json_validation_object = attribute.json_template[:validation]

      expect(json_validation_object.to_json).to eq("{\"type\":\"Regexp\",\"value\":\"(?-mix:^[a-zA-Z]{1}$)\"}")
    end

    it "includes attribute groups" do
      attribute = Magma::Attribute.new(
        attribute_name: "stench",
        attribute_group: "odors"
      )

      template = attribute.json_template

      expect(template[:attribute_group]).to eq("odors")
    end
  end

  describe "#revision_to_loader" do
    it "returns entry for editable string backed options" do
      attribute = Magma::StringAttribute.new(
        project_name: :project,
        attribute_name: "name",
        model_name: :model,
        description: "Old name"
      )

      entry = attribute.revision_to_loader(:description, "New name")

      expect(entry).to eq([:name, "New name"])
    end

    it "returns entry for editable JSON backed options" do
      attribute = Magma::MatchAttribute.new(
        project_name: :project,
        attribute_name: "name",
        model_name: :model,
        validation: { type: "Array", value: [1, 2, 3] }
      )

      entry = attribute.revision_to_loader(:validation, { type: "Array", value: [4, 5, 6] })

      expect(entry[0]).to eq(:name)
      expect(entry[1].to_json).to eq({type: "Array", value: [4, 5, 6]}.to_json)
    end
  end

  describe "#revision_to_payload" do
    it "returns entry for editable string backed options" do
      attribute = Magma::StringAttribute.new(
        project_name: :project,
        attribute_name: "name",
        model_name: :model,
        description: "Old name"
      )

      entry = attribute.revision_to_payload(
        :description,
        "New name",
        Etna::User.new({
          email: "outis@mountolympus.org"
        })
      )

      expect(entry).to eq([:name, "New name"])
    end

    it "returns entry for editable JSON backed options" do
      attribute = Magma::MatchAttribute.new(
        project_name: :project,
        attribute_name: "name",
        model_name: :model,
        validation: { type: "Array", value: [1, 2, 3] }
      )

      entry = attribute.revision_to_payload(
        :validation,
        { type: "Array", value: [4, 5, 6] },
        Etna::User.new({
          email: "outis@mountolympus.org"
        }))

      expect(entry[0]).to eq(:name)
      expect(entry[1].to_json).to eq({type: "Array", value: [4, 5, 6]}.to_json)
    end
  end

  describe "#query_to_payload" do
    it "returns attribute data value for a JSON payload" do
      attribute = Magma::Attribute.new(
        attribute_name: "name",
        description: "Old name")

      query = attribute.query_to_payload(
        Magma::Answer.new("New name")
      )

      expect(query).to eq("New name")
    end
  end

  describe "#query_to_tsv" do
    it "returns attribute value for a TSV payload" do
      attribute = Magma::Attribute.new(
        attribute_name: "name",
        description: "Old name")

      query = attribute.query_to_tsv("New name")

      expect(query).to eq("New name")
    end
  end

  describe "#validation_object" do
    it "builds ArrayValidationObjects using validation options from the database" do
      attribute = Magma::StringAttribute.new(
        project_name: "project",
        model_name: "model",
        attribute_name: "name",
        validation: { type: "Array", value: ["a", "b", "c"] }
      )

      expect(attribute.validation_object.validate("a")).to eq(true)
    end

    it "builds ArrayValidationObjects using validation options defined on the attribute" do
      attribute = Magma::Attribute.new(
        attribute_name: "name",
        validation: { type: "Array", value: ["a", "b", "c"] }
      )

      expect(attribute.validation_object.validate("a")).to eq(true)
    end

    it "builds RangeValidationObjects using validation options from the database" do
      attribute = Magma::IntegerAttribute.new(
        project_name: "project",
        model_name: "model",
        attribute_name: "name",
        validation: { type: "Range", begin: 1, end: 10 }
      )

      expect(attribute.validation_object.validate(5)).to eq(true)
    end

    it "builds RangeValidationObjects using validation options defined on the attribute" do
      attribute = Magma::Attribute.new(
        attribute_name: "name",
        validation: { type: "Range", begin: 1, end: 10 }
      )

      expect(attribute.validation_object.validate(5)).to eq(true)
    end

    it "builds RegexpValidationObjects using validation options from the database" do
      attribute = Magma::StringAttribute.new(
        project_name: "project",
        model_name: "model",
        attribute_name: "name",
        validation: { type: "Regexp", value: "^[a-zA-Z]{1}$" }
      )

      expect(attribute.validation_object.validate("A")).to eq(true)
    end

    it "builds RegexpValidationObjects using validation options defined on the attribute" do
      attribute = Magma::Attribute.new(
        attribute_name: "name",
        validation: { type: "Regexp", value: /^[a-zA-Z]{1}$/ }
      )

      expect(attribute.validation_object.validate("A")).to eq(true)
    end
  end

  describe "#valid?" do
    it "returns false if type doesn't match an existing Magma::Attribute" do
      attribute = Magma::Attribute.new(
        attribute_name: "name",
        type: "invalid"
      )

      expect(attribute.valid?).to eq(false)
    end
  end
end
