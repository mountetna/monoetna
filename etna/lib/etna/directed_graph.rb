class DirectedGraph
  def initialize
    @children = {}
    @parents = {}
  end

  attr_reader :children
  attr_reader :parents

  def include?(node)
    @children.include?(node)
  end

  def add_node(node)
    [@children[node] ||= {}, @parents[node] ||= {}]
  end

  def add_connection(parent, child)
    children, parent_parents = add_node(parent)
    parents, child_children = add_node(child)
    children[child] = child_children
    parents[parent] = parent_parents
  end

  def order_nodes
    processed = Set.new
    result = []

    queue = @children.keys.dup.sort.map { |k| [k, []] }
    while (node, parents = queue.shift)
      next if processed.include?(node)

      unseen = @children[node].keys.filter { |i| !processed.include?(i) }
      if unseen.empty?
        result << node
        processed.add(node)
        next
      end

      new_parents = parents.dup
      new_parents << node

      unseen.each do |unseen_node|
        if parents.include?(unseen_node)
          raise "Cannot order nodes, cyclical dependency found #{parents.join(' -> ')} -> #{unseen_node}"
        end

        queue.unshift([unseen_node, new_parents])
      end

      queue.unshift([node, parents])
    end

    result
  end

  def descendants(parent)
    seen = Set.new

    seen.add(parent)
    queue = @children[parent].keys.dup
    parent_queue = @parents[parent].keys.dup

    # Because this is not an acyclical graph, the definition of descendants needs to be stronger;
    # here we believe that any path that would move through --any-- parent to this child would not be considered
    # descendant, so we first find all those parents and mark them as 'seen' so that they are not traveled.
    while next_parent = parent_queue.pop
      next if seen.include?(next_parent)
      seen.add(next_parent)
      parent_queue.push(*@parents[next_parent].keys)
    end

    queue = queue.nil? ? [] : queue.dup
    paths = {}

    while child = queue.pop
      next if seen.include? child
      seen.add(child)
      path = (paths[child] ||= [parent])

      @children[child].keys.each do |child_child|
        queue.push child_child

        unless paths.include? child_child
          paths[child_child] = path + [child]
        end
      end
    end

    paths
  end
end