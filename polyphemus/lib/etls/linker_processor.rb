class Polyphemus::LinkerProcessor
  attr_reader :linker, :model_name

  def initialize(linker:, model_name:)
    @linker = linker
    @model_name = model_name
  end

  def process(cursor, files)
    if cursor[:project_name] != linker.project_name || cursor[:bucket_name] != linker.bucket_name
      raise "Invalid linker configuration: input cursor project and bucket do not match linker specification."
    end

    puts "linking #{files}"
    linker.link(
      model_name: model_name,
      files: files,
    )
  end
end