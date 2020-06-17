class BasicError < StandardError

  attr_reader :type
  attr_reader :id
  attr_reader :method

  def initialize(type, id, method)

    @type = type
    @id = id
    @method = method
  end
end