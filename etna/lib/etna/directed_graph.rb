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

    puts "adding connection #{child} < #{parent}"
  end

  def descendants(parent)
    seen = Set.new

    seen.add(parent)
    queue = @children[parent].keys.dup
    puts "initial keys #{queue}"
    parent_queue = @parents[parent].keys.dup

    while next_parent = parent_queue.pop
      next if seen.include?(next_parent)
      puts "ignoring #{next_parent} for exclusion as parent"
      seen.add(next_parent)
      parent_queue.push(*@parents[next_parent].keys)
    end

    queue = queue.nil? ? [] : queue.dup
    paths = {}

    while child = queue.pop
      next if seen.include? child
      puts "considering #{child}"
      seen.add(child)
      path = (paths[child] ||= [parent])

      @children[child].keys.each do |child_child|
        queue.push child_child

        unless paths.include? child_child
          paths[child_child] = path + [child]
        else
          puts "ignoring #{child_child}"
        end
      end
    end

    paths
  end
end