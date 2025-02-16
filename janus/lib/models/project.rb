class Project < Sequel::Model
  one_to_many :permissions

  PROJECT_NAME_MATCH=/(?!pg_)[a-z][a-z0-9]*(_[a-z][a-z0-9]*)*/

  def self.valid_name?(project_name)
    project_name =~ /\A#{PROJECT_NAME_MATCH.source}\Z/
  end

  PROJECT_TYPES={
    'team' => {
      public: false
    },
    'resource' => {
      public: true
    },
    'community' => {
      public: true,
      requires_agreement: true
    },
    'template' => {
      public: true
    }
  }

  PUBLIC_PROJECT_TYPES=PROJECT_TYPES.select do |type, props|
    props[:public]
  end.keys

  def requires_agreement?
    !!Project::PROJECT_TYPES[project_type][:requires_agreement]
  end

  def to_hash
    {
      project_name: project_name,
      project_name_full: project_name_full,
      project_type: project_type,
      permissions: permissions.map(&:to_hash),
      cc_text: cc_text,
      contact_email: contact_email
    }
  end
end
