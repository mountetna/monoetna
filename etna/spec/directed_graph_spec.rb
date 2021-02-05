describe DirectedGraph do
  describe '#paths_from' do
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
end