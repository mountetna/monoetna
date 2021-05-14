class DirectedGraph
  def initialize
    @children = {}
    @parents = {}
  end

  attr_reader :children
  attr_reader :parents

  def full_parentage(n)
    [].tap do |result|
      q = @parents[n].keys.dup
      seen = Set.new

      until q.empty?
        n = q.shift
        next if seen.include?(n)
        seen.add(n)

        result << n
        q.push(*@parents[n].keys)
      end

      result.uniq!
    end
  end

  def as_normalized_hash(root, include_root = true)
    q = [root]
    {}.tap do |result|
      if include_root
        result[root] = []
      end

      seen = Set.new

      until q.empty?
        n = q.shift
        next if seen.include?(n)
        seen.add(n)

        parentage = full_parentage(n)

        @children[n].keys.each do |child_node|
          q << child_node

          if result.include?(n)
            result[n] << child_node
          end

          parentage.each do |grandparent|
            result[grandparent] << child_node if result.include?(grandparent)
          end

          # Depending on the graph shape, diamonds could lead to
          #   resetting of previously calculated dependencies.
          # Here we avoid resetting existing entries in `result`
          #   and instead concatenate them if they already exist.
          result[child_node] = [] unless result.include?(child_node)
          result[n].concat(result[child_node]) if result.include?(child_node) && result.include?(n)
        end
      end

      result.values.each(&:uniq!)
    end
  end

  def add_connection(parent, child)
    children = @children[parent] ||= {}
    child_children = @children[child] ||= {}

    children[child] = child_children

    parents = @parents[child] ||= {}
    parent_parents = @parents[parent] ||= {}
    parents[parent] = parent_parents
  end

  def serialized_path_from(root, include_root = true)
    seen = Set.new
    [].tap do |result|
      result << root if include_root
      seen.add(root)
      path_q = paths_from(root, include_root)
      traversables = path_q.flatten

      until path_q.empty?
        next_path = path_q.shift
        next if next_path.nil?

        until next_path.empty?
          next_n = next_path.shift
          next if next_n.nil?
          next if seen.include?(next_n)

          if @parents[next_n].keys.any? { |p| !seen.include?(p) && traversables.include?(p) }
            next_path.unshift(next_n)
            path_q.push(next_path)
            break
          else
            result << next_n
            seen.add(next_n)
          end
        end
      end
    end
  end

  def paths_from(root, include_root = true)
    [].tap do |result|
      parents_of_map = descendants(root, include_root)
      seen = Set.new

      parents_of_map.to_a.sort_by { |k, parents| [-parents.length, k.inspect] }.each do |k, parents|
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

  def descendants(parent, include_root = true)
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
      path = (paths[child] ||= (include_root ? [parent] : []))

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