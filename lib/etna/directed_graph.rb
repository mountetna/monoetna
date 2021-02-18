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

  def paths_from(root)
    [].tap do |result|
      parents_of_map = descendants(root)
      seen = Set.new

      parents_of_map.to_a.sort_by { |k, parents| [-parents.length, k] }.each do |k, parents|
        unless seen.include?(k)
          if @children[k].keys.empty?
            result << parents + [k]
          else
            @children[k].keys.dup.sort.each do |c|
              result << parents + [k, c]
            end
          end
        end

        parents.each { |p| seen.add(p) }
      end
    end
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