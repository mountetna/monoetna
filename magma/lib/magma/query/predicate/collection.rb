class Magma
  class CollectionPredicate < Magma::Predicate
  # This object should just return the list of aggregated values
  #   for a collection.
    attr_reader :model

    def initialize question, model, *query_args
      super(question)
      @model = model
      require 'pry'
      binding.pry
      process_args(query_args)
    end

    verb '::all' do
      child :record_child

      extract do |table, return_identity|
        require 'pry'
        binding.pry
        table.map do |identifier,rows|
          next unless identifier
          [ identifier, rows ]
        end.compact
      end
    end

    verb '::first' do
      child :record_child

      extract do |table, return_identity|
        require 'pry'
        binding.pry
        table.group_by do |row|
          row[identity]
        end.first&.last
      end
    end

    def to_hash
      super.merge(
        model: model
      )
    end

    private

    def record_child
      RecordPredicate.new(@question, @model, alias_name, *@query_args)
    end
  end
end
