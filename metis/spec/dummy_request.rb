# A dummy class that behaves as a Rack::Request
class DummyRequest

  def initialize

    @req = {

      'directory'=> '/ipi/melanoma/',
      'expires'=> 600,
      'algorithm'=> 'MD5',
      'timestamp'=> 1473887352,
      'type'=> 'solid',
      'user_email'=> 'jason.cater@ucsf.edu',
      'auth_token'=> 'f976c1a092978d1562e6bb7b16d9a873',
      'file_name'=> 'rand.data.64B-ec65f526379d74e8e4dd60a56f37c868',
      'signature'=> '20e830f373664673b039239ccbfff5d2'
    }
  end

  def POST
    
    @req
  end

  def post?

    true
  end

  def munge

    @req['signature'] = ''
  end
end