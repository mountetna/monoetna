describe WorkflowError do
  describe 'mark_error! & find_error' do
    it 'works' do
      WorkflowError.mark_error!(
          hash: 'abcde', message: 'the message')

      expect(WorkflowError.find_error(hash: 'abcde')).to_not be_nil
      expect(WorkflowError.find_error(hash: 'abcde').message).to eql('the message')
      expect(WorkflowError.find_error(hash: 'abcdef')).to be_nil
    end

    it 'is idempotent' do
      expect(WorkflowError.find_error(hash: 'abcde')).to be_nil

      WorkflowError.mark_error!(
          hash: 'abcde', message: 'the message1')
      WorkflowError.mark_error!(
          hash: 'abcde', message: 'the message2')

      expect(WorkflowError.find_error(hash: 'abcde').message).to eql('the message1')
    end
  end

  describe 'clear!' do
    it 'is idempotent' do
      error = WorkflowError.mark_error!(hash: 'abcde', message: 'the message1')
      WorkflowError.mark_error!(hash: 'abcdef', message: 'other message')
      WorkflowError.mark_error!(hash: 'abcde', message: 'the message1')

      expect(WorkflowError.find_error(hash: 'abcde').message).to eql('the message1')
      expect(WorkflowError.find_error(hash: 'abcdef').message).to eql('other message')

      error.clear!
      expect(WorkflowError.find_error(hash: 'abcde')).to be_nil
      # Does not clear other errors.
      expect(WorkflowError.find_error(hash: 'abcdef').message).to eql('other message')

      # does not raise an exception to delete twice
      error.clear!
      error2 = WorkflowError.mark_error!(hash: 'abcde', message: 'the message2')

      # Can create new errors for the same hash.
      expect(WorkflowError.find_error(hash: 'abcde').message).to eql('the message2')

      # Does not clear errors by a different uuid
      error.clear!
      expect(WorkflowError.find_error(hash: 'abcde').message).to eql('the message2')

      error2.clear!
      expect(WorkflowError.find_error(hash: 'abcde')).to be_nil
    end
  end
end