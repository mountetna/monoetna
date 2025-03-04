describe WithSlackNotifications do
  class TestNotifier
    include WithSlackNotifications
  end

  let(:url) {Polyphemus.instance.config(:slack_webhook_url)}

  it 'sends a mesage to slack' do
    stub_slack
    msg = 'Hail the storm-footed, golden-winged daughter of Thaumas'

    TestNotifier.new.notify_slack( msg, channel: '#agora')
    expect(WebMock).to have_requested(:post, url).with(
      body: hash_including(text: msg)
    )
  end
end
