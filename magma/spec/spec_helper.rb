require 'yaml'
require 'logger'
require 'factory_bot'
require 'database_cleaner'
require 'rack/test'
require 'simplecov'
require 'timecop'
require 'webmock/rspec'

SimpleCov.start
require 'bundler'
Bundler.require(:default, :test)

require 'etna'
require 'etna/spec/event_log'


ENV['MAGMA_ENV'] = 'test'

# !!!ACHTUNG!!!
# We have to build the application from scratch, since config.ru may be
# wrong. Check to ensure that this is not out-of-date in some important way!

require_relative '../lib/magma/server'
require_relative '../lib/magma'

Magma.instance.configure(YAML.load(File.read('config.yml')))

# This is a bit of a hack because we don't yet fully support all attribute features in the csv, but the csv is also the
# most reliable way to force migrations and project setup.  Furthermore, sequel requires object full reloads when performing
# certain kinds of reflection + migrations during runtime, so we're forced to flush and reload a bunch all over the place.
# This can get nicer with deeper plumbing of the way sequel is used in magma and improvements to the csv interface.
def reset_labors_project
  Magma.instance.db[:attributes].truncate
  Magma.instance.db[:models].truncate

  # Remove existing labors schema, start scratch every time.
  Magma.instance.db.run <<-SQL
DO $$
BEGIN
    IF EXISTS(
        SELECT schema_name
          FROM information_schema.schemata
          WHERE schema_name = 'labors'
      )
    THEN
      EXECUTE 'DROP SCHEMA labors CASCADE';
    END IF;
END
$$;
  SQL
end

def load_labors_project
  stub_event_log
  reset_labors_project
  WebMock.enable!

  Magma.instance.setup_sequel
  Magma.instance.setup_logger
  Magma.instance.magma_projects.clear
  Magma.instance.load_db_projects

  magma_client = Etna::Clients::LocalMagmaClient.new

  magma_client.update_model(Etna::Clients::Magma::UpdateModelRequest.new(
    project_name: "labors",
    actions: [Etna::Clients::Magma::AddProjectAction.new(no_metis_bucket: true)]))

  # This initializes and migrates the schema.
  workflow = Etna::Clients::Magma::AddProjectModelsWorkflow.new(magma_client: magma_client)
  errors = []

  changeset = File.open("./spec/fixtures/labors_models_project_tree.csv", 'r') do |f|
    workflow.prepare_changeset_from_csv(io: f) do |err|
      errors << err
    end
  end

  raise "Failed to load labors project, found errors in input csv: #{errors.join('\n')}" unless errors.empty?
  sync_workflow = workflow.plan_synchronization(changeset, "labors")
  print "Adding models and attributes"
  sync_workflow.update_block = Proc.new do |action|
    print "."
    Object.class_eval { remove_const(:Labors) if Object.const_defined?(:Labors) }
    Magma.instance.magma_projects.clear
    Magma.instance.load_models(false)
  end
  puts
  sync_workflow.execute_planned!

  Object.class_eval { remove_const(:Labors) if Object.const_defined?(:Labors) }
  Magma.instance.magma_projects.clear
  Magma.instance.load_models(false)

  # This side loads some attributes, specifically the dictionary attribute, which is not currently supported
  # via the csv api.
  YAML.load(File.read("./spec/fixtures/labors_model_attributes.yml")).each do |model_name, attributes|
    Magma.instance.db[:models].where(
      project_name: "labors",
      model_name: model_name,
    ).update(
      dictionary: attributes.delete("dictionary")
    )

    attributes.each do |attribute_name, options|
      row = options
      row["column_name"] = attribute_name unless row["column_name"]

      Magma.instance.db[:attributes].where(
        project_name: "labors",
        model_name: model_name,
        attribute_name: attribute_name,
      ).update(row)
    end
  end

  Object.class_eval { remove_const(:Labors) if Object.const_defined?(:Labors) }
  Magma.instance.magma_projects.clear
  Magma.instance.load_models(false)

  require_relative './labors/metrics/labor_metrics'
end

Magma.instance.setup_db
load_labors_project

OUTER_APP = Rack::Builder.new do
  use Etna::ParseBody
  use Etna::SymbolizeParams

  use Etna::TestAuth

  run Magma::Server.new
end

# This file was generated by the `rspec --init` command. Conventionally, all
# specs live under a `spec` directory, which RSpec adds to the `$LOAD_PATH`.
# The generated `.rspec` file contains `--require spec_helper` which will cause
# this file to always be loaded, without a need to explicitly require it in any
# files.
#
# Given that it is always loaded, you are encouraged to keep this file as
# light-weight as possible. Requiring heavyweight dependencies from this file
# will add to the boot time of your test suite on EVERY test run, even for an
# individual file that may not need all of that loaded. Instead, consider making
# a separate helper file that requires the additional dependencies and performs
# the additional setup, and require it from the spec files that actually need
# it.
#
# See http://rubydoc.info/gems/rspec-core/RSpec/Core/Configuration
RSpec.configure do |config|
  # rspec-expectations config goes here. You can use an alternate
  # assertion/expectation library such as wrong or the stdlib/minitest
  # assertions if you prefer.
  config.expect_with :rspec do |expectations|
    expectations.on_potential_false_positives = :nothing
    # This option will default to `true` in RSpec 4. It makes the `description`
    # and `failure_message` of custom matchers include text for helper methods
    # defined using `chain`, e.g.:
    #     be_bigger_than(2).and_smaller_than(4).description
    #     # => "be bigger than 2 and smaller than 4"
    # ...rather than:
    #     # => "be bigger than 2"
    expectations.include_chain_clauses_in_custom_matcher_descriptions = true
  end

  # suppress annoying console output during the run
  original_stderr = $stderr
  original_stdout = $stdout
  config.before(:all) do
    #$stderr = File.open(File::NULL, "w")
    #$stdout = File.open(File::NULL, "w")
  end
  config.after(:all) do
    $stderr = original_stderr
    $stdout = original_stdout
  end

  # rspec-mocks config goes here. You can use an alternate test double
  # library (such as bogus or mocha) by changing the `mock_with` option here.
  config.mock_with :rspec do |mocks|
    # Prevents you from mocking or stubbing a method that does not exist on
    # a real object. This is generally recommended, and will default to
    # `true` in RSpec 4.
    mocks.verify_partial_doubles = true
  end

  # This option will default to `:apply_to_host_groups` in RSpec 4 (and will
  # have no way to turn it off -- the option exists only for backwards
  # compatibility in RSpec 3). It causes shared context metadata to be
  # inherited by the metadata hash of host groups and examples, rather than
  # triggering implicit auto-inclusion in groups with matching metadata.
  config.shared_context_metadata_behavior = :apply_to_host_groups

# The settings below are suggested to provide a good initial experience
# with RSpec, but feel free to customize to your heart's content.
=begin
  # This allows you to limit a spec run to individual examples or groups
  # you care about by tagging them with `:focus` metadata. When nothing
  # is tagged with `:focus`, all examples get run. RSpec also provides
  # aliases for `it`, `describe`, and `context` that include `:focus`
  # metadata: `fit`, `fdescribe` and `fcontext`, respectively.
  config.filter_run_when_matching :focus

  # Allows RSpec to persist some state between runs in order to support
  # the `--only-failures` and `--next-failure` CLI options. We recommend
  # you configure your source control system to ignore this file.
=end
  config.example_status_persistence_file_path = "spec/examples.txt"
=begin
  # Limits the available syntax to the non-monkey patched syntax that is
  # recommended. For more details, see:
  #   - http://rspec.info/blog/2012/06/rspecs-new-expectation-syntax/
  #   - http://www.teaisaweso.me/blog/2013/05/27/rspecs-new-message-expectation-syntax/
  #   - http://rspec.info/blog/2014/05/notable-changes-in-rspec-3/#zero-monkey-patching-mode
  config.disable_monkey_patching!

  # This setting enables warnings. It's recommended, but in some cases may
  # be too noisy due to issues in dependencies.
  config.warnings = true

  # Many RSpec users commonly either run the entire suite or an individual
  # file, and it's useful to allow more verbose output when running an
  # individual spec file.
  if config.files_to_run.one?
    # Use the documentation formatter for detailed output,
    # unless a formatter has already been configured
    # (e.g. via a command-line flag).
    config.default_formatter = "doc"
  end

  # Print the 10 slowest examples and example groups at the
  # end of the spec run, to help surface which specs are running
  # particularly slow.
  config.profile_examples = 10

  # Run specs in random order to surface order dependencies. If you find an
  # order dependency and want to debug it, you can fix the order by providing
  # the seed, which is printed after each run.
  #     --seed 1234

  # Seed global randomization in this process using the `--seed` CLI option.
  # Setting this allows you to use `--seed` to deterministically reproduce
  # test failures related to randomization by passing the same `--seed` value
  # as the one that triggered the failure.
  Kernel.srand config.seed
=end

  config.order = :random

  config.include FactoryBot::Syntax::Methods

  config.before(:suite) do
    FactoryBot.find_definitions
    #DatabaseCleaner[:sequel].db = Magma.instance.db
    #DatabaseCleaner.strategy = :transaction
    #DatabaseCleaner.clean_with(:truncation, except: ["models", "attributes", "schema_info"])
  end

  config.around(:each) do |example|
    #DatabaseCleaner.cleaning { example.run }
    Magma.instance.db.transaction(:rollback=>:always, :auto_savepoint=>true){example.run}
  end
end

FactoryBot.define do
  factory :labor, class: Labors::Labor do
    to_create(&:save)
    sequence(:name) { |n| "labor#{n}" }
    sequence(:number) { |n| n+1 }
    sequence(:completed) { |n| [2, 5].include?(number) ? false : true }

    trait :lion do
      name { 'Nemean Lion' }
      number { 1 }
      completed { true }
    end

    trait :hydra do
      name { 'Lernean Hydra' }
      number { 2 }
      completed { false }
    end

    trait :stables do
      name { 'Augean Stables' }
      number { 5 }
      completed { false }
    end

    trait :hind do
      name { 'Ceryneian Hind' }
      number { 3 }
      completed { true }
    end
  end

  factory :monster, class: Labors::Monster do
    to_create(&:save)
    sequence(:name) { |n| "monster#{n}" }

    trait :lion do
      name { 'Nemean Lion' }
      species { 'lion' }
    end

    trait :hydra do
      name { 'Lernean Hydra' }
      species { 'hydra' }
    end

    trait :hind do
      name { 'Ceryneian Hind' }
      species { 'red deer' }
    end

    trait :birds do
      name { 'Stymphalian Birds' }
      species { 'marsh bird' }
    end
  end

  factory :victim, class: Labors::Victim do
    to_create(&:save)
    sequence(:name) { |n| "victim#{n}" }
  end

  factory :sidekick, class: Labors::Sidekick do
    to_create(&:save)
    sequence(:name) { |n| "sidekick#{n}" }
  end

  factory :prize, class: Labors::Prize do
    to_create(&:save)
    sequence(:name) { |n| "prize#{n}" }
  end

  factory :codex, class: Labors::Codex do
    to_create(&:save)

    trait :lion do
      monster { 'Nemean Lion' }
    end

    trait :hydra do
      monster { 'Lernean Hydra' }
    end
  end
  factory :aspect, class: Labors::Aspect do
    to_create(&:save)
  end

  factory :project, class: Labors::Project do
    to_create(&:save)
  end

  factory :habitat, class: Labors::Habitat do
    to_create(&:save)
  end

  factory :characteristic, class: Labors::Characteristic do
    to_create(&:save)
    sequence(:name) { |n| "characteristic#{n}" }
  end

  factory :wound, class: Labors::Wound do
    to_create(&:save)
  end

  factory :grammar, class: Magma::Gnomon::Grammar do
    to_create(&:save)
  end

  factory :identifier, class: Magma::Gnomon::Identifier do
    to_create(&:save)
  end

  factory :flag, class: Magma::Flag do
    to_create(&:save)

    trait :gnomon_none do
      project_name {"labors"}
      flag_name { Magma::Flags::GNOMON_MODE[:name]}
      value {Magma::Flags::GNOMON_MODE[:none]}
    end

    trait :gnomon_identifier do
      project_name {"labors"}
      flag_name { Magma::Flags::GNOMON_MODE[:name]}
      value {Magma::Flags::GNOMON_MODE[:identifier]}
    end

    trait :gnomon_pattern do
      project_name {"labors"}
      flag_name { Magma::Flags::GNOMON_MODE[:name]}
      value {Magma::Flags::GNOMON_MODE[:pattern]}
    end

  end

end

def fixture(name)
  File.join(File.dirname(__FILE__), "fixtures/#{name}.txt")
end

def json_body(response=nil)
  JSON.parse((response || last_response).body, symbolize_names: true)
end

def json_document model, record_name
  json_body[:models][model.to_sym][:documents][record_name.to_sym]
end

AUTH_USERS = {
  superuser: {
    email: 'zeus@twelve-labors.org', name: 'Zeus', perm: 'a:administration'
  },
  admin: {
    email: 'hera@twelve-labors.org', name: 'Hera', perm: 'a:labors'
  },
  editor: {
    email: 'eurystheus@twelve-labors.org', name: 'Eurystheus', perm: 'e:labors'
  },
  privileged_editor: {
    email: 'copreus@twelve-labors.org', name: 'Copreus', perm: 'E:labors'
  },
  viewer: {
    email: 'hercules@twelve-labors.org', name: 'Hercules', perm: 'v:labors'
  },
  non_user: {
    email: 'nessus@centaurs.org', name: 'Nessus', perm: ''
  },
  guest: {
    email: 'sinon@troy.org', name: 'Sinon', perm: 'g:labors'
  }
}

def auth_header(user_type)
  header(*Etna::TestAuth.token_header(AUTH_USERS[user_type]))
end

def json_post(endpoint, hash)
  post(
    "/#{endpoint.to_s.reverse.chomp('/').reverse}",
    hash.to_json,
    {'CONTENT_TYPE'=> 'application/json'}
  )
end

def setup_metis_bucket_stubs(project_name)
  route_payload = JSON.generate([
    {:method=>"POST", :route=>"/:project_name/bucket/create/:bucket_name", :name=>"bucket_create", :params=>["project_name","bucket_name"]}
  ])
  stub_request(:options, 'https://metis.test').
  to_return(status: 200, body: route_payload, headers: {'Content-Type': 'application/json'})

  route_payload = JSON.generate([
    {:bucket=>{}}
  ])
  stub_request(:post, /https:\/\/metis.test\/#{project_name}\/bucket\/create/).
    to_return(status: 200, body: route_payload, headers: {'Content-Type': 'application/json'})
end

def stub_validation(model, att_name, new_validation)
  validation_stubs[model] ||= {}
  validation_stubs[model][att_name] = model.attributes[att_name].validation
  model.attributes[att_name].validation = new_validation
end

def remove_validation_stubs
  validation_stubs.each do |model,atts|
    atts.each do |att_name, old_validation|
      model.attributes[att_name].validation = old_validation
    end
  end
end

def set_date_shift_root(model_name, value)
  Magma::SetDateShiftRootAction.new("labors", {
    action_name: "set_date_shift_root",
    model_name: model_name,
    date_shift_root: value,
  }).perform
end

def stub_date_shift_data(project)
  hydra = create(:labor, :hydra, project: project)
  lion = create(:labor, :lion, project: project)
  @hind = create(:labor, :hind, project: project)

  @lion_monster = create(:monster, :lion, labor: lion)
  @hydra_monster = create(:monster, :hydra, labor: hydra)

  @john_doe = create(:victim, name: 'John Doe', monster: @lion_monster, country: 'Italy')
  jane_doe = create(:victim, name: 'Jane Doe', monster: @lion_monster, country: 'Greece')

  @susan_doe = create(:victim, name: 'Susan Doe', monster: @hydra_monster, country: 'Italy')
  shawn_doe = create(:victim, name: 'Shawn Doe', monster: @hydra_monster, country: 'Greece')

  @john_arm = create(:wound, victim: @john_doe, location: 'Arm', severity: 5)
  create(:wound, victim: @john_doe, location: 'Leg', severity: 1)
  create(:wound, victim: jane_doe, location: 'Arm', severity: 2)
  create(:wound, victim: jane_doe, location: 'Leg', severity: 4)
  @susan_arm = create(:wound, victim: @susan_doe, location: 'Arm', severity: 3)
  create(:wound, victim: @susan_doe, location: 'Leg', severity: 3)
  create(:wound, victim: shawn_doe, location: 'Arm', severity: 1)
  create(:wound, victim: shawn_doe, location: 'Leg', severity: 1)
end

def iso_date_str(value)
  DateTime.parse(value).iso8601
end

# Gnomon helpers
def create_identifier(id, params={})
  identifier = create(
    :identifier, {
    project_name: 'labors',
    author: "Hera|hera@twelve-labors.org",
    identifier: id,
  }.merge(params)
  )
end

def register_flag(flag_hash)
    # This dynamically creates constants in the Flags mode and "registers" them
    # The first argument is just the name of the constant, it can be random.
    constant_name = 'TEST_' + (1..4).map { ('A'..'Z').to_a.sample }.join
    Magma::Flags.const_set(constant_name, flag_hash)
end

def unregister_flags
  Magma::Flags.constants.select { |const| const.to_s.start_with?('TEST_') }.each do |const|
    Magma::Flags.module_eval { remove_const(const) }
  end
end

def create_flags_in_db(flag_hash)
  flag_hash.each do |k,v|
    create(:flag, project_name: "labors", flag_name: k, value: v)
  end
end

VALID_GRAMMAR_CONFIG={
  tokens: {
    PROJECT: {
      label: "project",
      values: {
        "The Twelve Labors of Hercules": "The Twelve Labors of Hercules"
      }
    },
    PROJ: {
      label: "project",
      values: {
        "LABORS": "The Twelve Labors of Hercules"
      }
    },
    LABOR: {
      label: "labor",
      values: {
        "The Nemean Lion": "The Nemean Lion",
        "The Lernean Hydra": "The Lernean Hydra"
      }
    },
    LAB: {
      label: "labor",
      values: {
        "LION": "The Nemean Lion",
        "HYDRA": "The Lernean Hydra"
      }
    },
    VILL: {
      label: "Village type",
      values: {
        "V": "Village",
        "H": "Hamlet"
      }
    },
    VICT: {
      label: "Victim type",
      values: {
        "S": "Soldier",
        "C": "Civilian"
      }
    },
    SEP: {
      label: "Separator",
      values: {
        "-": "# Separator"
      }
    }

  },
  synonyms: [
    [ "PROJ", "PROJECT" ],
    [ "LAB", "LABOR" ]
  ],
  rules: {
    project: "PROJECT",
    labor: "LABOR",
    village: "PROJ SEP LAB SEP VILL .n",
    victim: ".village SEP VICT .n"
  }
}

HIERARCHY_GRAMMAR_CONFIG={
  tokens: {
    PROJECT: {
      label: "project",
      values: {
        "The Twelve Labors of Hercules": "The Twelve Labors of Hercules"
      }
    },
    PROJ: {
      label: "project",
      values: {
        "LABORS": "The Twelve Labors of Hercules"
      }
    },
    LABOR: {
      label: "labor",
      values: {
        "The Nemean Lion": "The Nemean Lion",
        "The Lernean Hydra": "The Lernean Hydra"
      }
    },
    LAB: {
      label: "labor",
      values: {
        "LION": "The Nemean Lion",
        "HYDRA": "The Lernean Hydra"
      }
    },
    VILL: {
      label: "Village type",
      values: {
        "V": "Village",
        "H": "Hamlet"
      }
    },
    MONSTER: {
      label: "monster",
      values: {
        "Nemean Lion": "Nemean Lion"
      }
    },
    MONST: {
      label: "monster",
      values: {
        "NEMEAN": "Nemean Lion"
      }
    },
    VICT: {
      label: "Victim type",
      values: {
        "S": "Soldier",
        "C": "Civilian"
      }
    },
    SEP: {
      label: "Separator",
      values: {
        "-": "# Separator"
      }
    }

  },
  synonyms: [
    [ "PROJ", "PROJECT" ],
    [ "LAB", "LABOR" ],
    [ "MONST", "MONSTER" ]
  ],
  rules: {
    project: "PROJECT",
    labor: "LABOR",
    monster: "PROJ SEP LAB SEP MONST",
    village: "PROJ SEP LAB SEP MONST SEP VILL .n",
    victim: ".village SEP VICT .n"
  }
}
