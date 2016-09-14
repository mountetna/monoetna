# General Utilities for Metis

module Utils

  # Rebuilds the request object from POST parameters and orders the values in 
  # an array for hashing.
  def Utils.generate_request(params)

    Conf::SIGNATURE_ITEMS.map do |item|
      
      params[item]
    end
  end

  # Takes an ordered array of request values and returns a signed hash.
  def Utils.sign_request(request, algo)

    signature = case algo.downcase
    when 'md5'    then sign_with_MD5(request)
    when 'sha256' then sign_with_SHA256(request)
    else ''
    end
  end

  # Takes an ordered array of request values, strigifies it, concatenates a
  # secret and hashes the resultant string with MD5.
  def Utils.sign_with_MD5(request)

    request_str = stringify_request(request) + Conf::SECRET_KEY

    md5 = Digest::MD5.new
    md5.update request_str
    md5.hexdigest
  end

  # Takes an ordered array of request values, strigifies it, concatenates a
  # secret and hashes the resultant string with SHA256.
  def Utils.sign_with_SHA256(request)

    request_str = stringify_request(request) + Conf::SECRET_KEY

    sha256 = Digest::SHA256.new
    sha256.update request_str
    sha256.hexdigest
  end

  # Takes an ordered array of request vaules and strigifies it. 
  def Utils.stringify_request(ordered_request)

    ordered_request.join
  end
end