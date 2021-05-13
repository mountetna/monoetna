describe DirectedGraph do
  let(:graph_one) do
    graph = DirectedGraph.new
    graph.add_connection('a', 'b')
    graph.add_connection('a', 'c')
    graph.add_connection('b', 'd')
    graph.add_connection('b', 'e')
    graph.add_connection('c', 'f')
    graph.add_connection('e', 'g')
    graph
  end

  let(:graph_two) do
    graph = DirectedGraph.new
    graph.add_connection('a', 'b')
    graph.add_connection('a', 'c')
    graph.add_connection('c', 'd')
    graph.add_connection('b', 'd')
    graph.add_connection('d', 'e')
    graph.add_connection('d', 'f')
    graph.add_connection('e', 'g')
    graph.add_connection('f', 'g')
    graph.add_connection('g', 'h')
    graph.add_connection('g', 'i')
    graph.add_connection('t', 'i')
    graph
  end

  let(:graph_three) do
    graph = DirectedGraph.new
    graph.add_connection('a', 'b')
    graph.add_connection('a', 'c')
    graph.add_connection('c', 'd')
    graph.add_connection('d', 'b')
    graph.add_connection('b', 'e')
    graph.add_connection('a', 'f')
    graph.add_connection('e', 'f')
    graph
  end

  describe '#paths_from' do
    describe 'divergent paths' do
      it 'works' do
        expect(graph_one.paths_from('a')).to eql([
            ["a", "b", "e", "g"],
            ["a", "b", "d"],
            ["a", "c", "f"],
        ])

        expect(graph_one.paths_from('a', false)).to eql([
            ["b", "e", "g"],
            ["b", "d"],
            ["c", "f"],
        ])
      end
    end

    describe 'with convergent and divergent paths' do
      it 'works' do
        expect(graph_two.paths_from('a')).to eql([
            ["a", "c", "d", "f", "g", "h"],
            ["a", "c", "d", "f", "g", "i"],
            ["a", "c", "d", "e", "g"],
            ["a", "b", "d"]
        ])
      end
    end
  end

  describe '#serialized_path_from' do
    it 'generates a linear serialization of convergent and divergent paths' do
      expect(graph_two.serialized_path_from('a')).to eql([
          "a", "c", "b", "d", "f", "e", "g", "h", "i"
      ])

      expect(graph_two.serialized_path_from('a', false)).to eql([
          "c", "b", "d", "f", "e", "g", "h", "i"
      ])
    end

    describe '#as_normalized_hash' do
      it 'works' do
        expect(graph_two.as_normalized_hash('a')).to eql({
            "a" => ["b", "c", "d", "e", "f", "g", "h", "i"],
            "b" => ["d", "e", "f", "g", "h", "i"],
            "c" => ["d", "e", "f", "g", "h", "i"],
            "d" => ["e", "f", "g", "h", "i"],
            "e" => ["g", "h", "i"],
            "f" => ["g", "h", "i"],
            "g" => ["h", "i"],
            "h" => [],
            "i" => [],
        })
        expect(graph_two.as_normalized_hash('a', false)).to eql({
            "b" => ["d", "e", "f", "g", "h", "i"],
            "c" => ["d", "e", "f", "g", "h", "i"],
            "d" => ["e", "f", "g", "h", "i"],
            "e" => ["g", "h", "i"],
            "f" => ["g", "h", "i"],
            "g" => ["h", "i"],
            "h" => [],
            "i" => [],
        })
        expect(graph_three.as_normalized_hash('a')).to eql({
          "a" => ["b", "c", "f", "e", "d"],
          "b" => ["e", "f"],
          "c" => ["e", "d", "f", "b"],
          "d" => ["f", "b", "e"],
          "e" => ["f"],
          "f" => [],
      })
      end
    end
  end
end