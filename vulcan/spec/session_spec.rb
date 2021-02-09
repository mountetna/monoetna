describe Session do
  describe '#as_json <-> from_json' do
    it 'works' do
      session = Session.new('project_name', 'workflow', 'thekey', {"a" => 1, "b" => 2})
      expect(Session.from_json(JSON.parse(session.as_json.to_json)).as_json).to eql(session.as_json)
    end
  end
end
