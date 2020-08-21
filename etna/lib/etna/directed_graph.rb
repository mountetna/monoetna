class DirectedGraph
  def initialize
    @children = {}
    @parents = {}
  end

  attr_reader :children
  attr_reader :parents

  def add_connection(parent, child)
    children = @children[parent] ||= {}
    child_children = @children[child] ||= {}

    children[child] = child_children

    parents = @parents[child] ||= {}
    parent_parents = @parents[parent] ||= {}
    parents[parent] = parent_parents
  end

  def descendants(parent)
    result = []
    seen = Set.new
    queue = @children[parent].keys
    queue = queue.nil? ? [] : queue.dup
    while child = queue.pop
      next if seen.include? child
      seen.add(child)
      result << child
      queue.push *@children[child].keys
    end

    result
  end
end