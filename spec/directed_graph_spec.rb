describe DirectedGraph do
  describe '#paths_from' do
    describe 'divergent paths' do
      it 'works' do
        graph = DirectedGraph.new
        graph.add_connection('a', 'b')
        graph.add_connection('a', 'c')
        graph.add_connection('b', 'd')
        graph.add_connection('b', 'e')
        graph.add_connection('c', 'f')
        graph.add_connection('e', 'g')

        expect(graph.paths_from('a')).to eql([
            ["a", "b", "e", "g"],
            ["a", "b", "d"],
            ["a", "c", "f"],
        ])
      end
    end

    describe 'with convergent and divergent paths' do
      it 'works' do
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

        expect(graph.paths_from('a')).to eql([
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

      expect(graph.serialized_path_from('a')).to eql([
          "a", "c", "b", "d", "f", "e", "g", "h", "i"
      ])
    end
  end
end