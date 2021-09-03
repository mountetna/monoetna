require_relative "../metis_folder_etl"

class Polyphemus::MetisFolderFilteringBaseEtl < Polyphemus::MetisFolderEtl
  def initialize(project_bucket_pairs:, folder_path_regexes: {}, limit: 20)
    @folder_path_regexes = folder_path_regexes
    super(
      project_bucket_pairs: project_bucket_pairs,
      limit: limit,
    )
  end

  def process(cursor, folders)
    raise "Subclass should do its filtering here"
  end

  def project_bucket_symbol(cursor)
    "#{cursor[:project_name]}_#{cursor[:bucket_name]}".to_sym
  end

  def filter_target_folders(cursor, folders)
    folders.select do |folder|
      folder_path_satisfies_regex?(cursor, folder)
    end
  end

  def should_skip_folder?(folder)
    # Subclasses can override
    false
  end

  def folder_path_satisfies_regex?(cursor, folder)
    project_bucket = project_bucket_symbol(cursor)
    regex = @folder_path_regexes[project_bucket]

    raise Polyphemus::EtlError.new("No regex for project / bucket: #{project_bucket}") if regex.nil?

    folder.folder_path =~ Regexp.new(regex)
  end
end
