require_relative "../lib/etna/permissions"

describe Etna::Permissions do
  it "can create from encoded permissions" do
    encoded_permissions = "a:labors;e:argo,olympics;v:constellations"
    perms = Etna::Permissions.from_encoded_permissions(encoded_permissions)

    expect(perms.to_string).to eq(encoded_permissions)
  end

  it "can create from hash" do
    permissions_hash = {
      "labors" => {
        role: :admin,
        restricted: false,
      },
      "olympics" => {
        role: :editor,
        restricted: true,
      },
      "argo" => {
        role: :editor,
        restricted: true,
      },
      "constellations" => {
        role: :viewer,
        restricted: false,
      },
    }
    perms = Etna::Permissions.from_hash(permissions_hash)

    expect(perms.to_hash).to eq(permissions_hash)
  end

  context "class instances" do
    before(:each) do
      @perms = Etna::Permissions.new([
        Etna::Permission.new("a", "labors"),
        Etna::Permission.new("E", "olympics"),
        Etna::Permission.new("E", "argo"),
        Etna::Permission.new("v", "constellations"),
      ])
    end

    it "serializes to string" do
      expect(@perms.to_string).to eq("E:argo,olympics;a:labors;v:constellations")
    end

    it "can convert to hash" do
      expect(@perms.to_hash).to eq({
        "labors" => {
          role: :admin,
          restricted: false,
        },
        "olympics" => {
          role: :editor,
          restricted: true,
        },
        "argo" => {
          role: :editor,
          restricted: true,
        },
        "constellations" => {
          role: :viewer,
          restricted: false,
        },
      })
    end

    it "can add a new permission" do
      @perms.add_permission(
        Etna::Permission.new("V", "fire")
      )

      expect(@perms.to_string).to eq("E:argo,olympics;V:fire;a:labors;v:constellations")
    end

    it "does not add permission if one exists for the given project" do
      @perms.add_permission(
        Etna::Permission.new("v", "argo")
      )

      expect(@perms.to_string).to eq("E:argo,olympics;a:labors;v:constellations")
    end
  end
end
