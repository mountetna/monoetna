class Magma
  module WithFilePredicateMethods
    attr_reader :requested_file_paths

    def initialize question, model, alias_name, attribute, *query_args
      super
      @md5_set = Md5Set.new(@question.user, @model)
      @updated_at_set = UpdatedAtSet.new(@question.user, @model)
    end
  end
end