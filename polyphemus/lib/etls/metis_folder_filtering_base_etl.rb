require_relative "../metis_folder_etl"

class Polyphemus::MetisFolderFilteringBaseEtl < Polyphemus::MetisFolderEtl
  def initialize(project_bucket_pairs:, folder_path_regexes: {}, limit: 20)
    @folder_path_regexes = folder_path_regexes
    super(
      project_bucket_pairs: project_bucket_pairs,
      limit: limit,
    )
  end

  def project_bucket_symbol(project_name:, bucket_name:)
    "#{project_name}_#{bucket_name}".to_sym
  end

  def filter_target_folders(cursor, folders)
    folders = folders.select do |folder|
      folder_path_satisfies_regex?(cursor, folder)
    end

    folders
  end

  def folder_path_satisfies_regex?(cursor, folder)
    return true unless cursor[:project_name] && cursor[:bucket_name]
    project_bucket = project_bucket_symbol(project_name: cursor[:project_name], bucket_name: cursor[:bucket_name])
    regex = @folder_path_regexes[project_bucket]

    raise Polyphemus::EtlError.new("No regex for project / bucket: #{project_bucket}") if regex.nil?

    folder.folder_path =~ Regexp.new(regex)
  end
end
