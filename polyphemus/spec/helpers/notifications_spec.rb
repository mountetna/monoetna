describe WithSlackNotifications do
  class TestNotifier
    include WithSlackNotifications
  end

  let(:url) {Polyphemus.instance.config(:slack_webhook_url)}
  let(:msg) {
    "Hail the storm-footed, golden-winged daughter of Thaumas\n" +
    "Shining in the raiment of her radiant spectrum"
  }

  it 'sends a mesage to slack' do
    stub_slack

    TestNotifier.new.notify_slack( msg, channel: '#agora')
    expect(WebMock).to have_requested(:post, url).with(
      body: hash_including(text: msg)
    )
  end

  context "message size" do
    before do
      WithSlackNotifications.send(:remove_const,:MESSAGE_SIZE)
      WithSlackNotifications.const_set(:MESSAGE_SIZE, 80)
    end

    after do
      WithSlackNotifications.send(:remove_const,:MESSAGE_SIZE)
      WithSlackNotifications.const_set(:MESSAGE_SIZE, 8000)
    end

    it 'splits a large message into several' do
      stub_slack

      TestNotifier.new.notify_slack( msg, channel: '#agora')
      expect(WebMock).to have_requested(:post, url).times(2)
    end
  end
end
