# sign_service.rb
# General signing utilities.

module SignService

  # Rebuilds the request object from POST parameters and orders the values in 
  # an array for hashing.
  def SignService.order_params(params)

    Conf::SIGNATURE_ITEMS.map do |item|

      params[item]
    end
  end

  # Verify that the appropriate parameters exsit for hashing or general usage.
  def SignService.request_parameters_valid?(params)

    valid = true

    Conf::SIGNATURE_ITEMS.each do |item|
      
      if !params.key?(item)

        valid = false
      end
    end

    return valid
  end

  # Takes an ordered array of request values and returns a signed hash.
  def SignService.sign_request(params, algo)

    signature = case algo.downcase
    when 'md5'    then sign_with_MD5(params)
    when 'sha256' then sign_with_SHA256(params)
    else ''
    end
  end

  # Takes an ordered array of request values, strigifies it, concatenates a
  # secret and hashes the resultant string with MD5.
  def SignService.sign_with_MD5(params)

    param_str = stringify_params(params) + Secrets::SECRET_KEY
    md5 = Digest::MD5.new
    md5.update param_str
    md5.hexdigest
  end

  # Takes an ordered array of request values, strigifies it, concatenates a
  # secret and hashes the resultant string with SHA256.
  def SignService.sign_with_SHA256(params)

    param_str = stringify_params(params) + Secrets::SECRET_KEY
    sha256 = Digest::SHA256.new
    sha256.update param_str
    sha256.hexdigest
  end

  # Takes an ordered array of request vaules and strigifies it. 
  def SignService.stringify_params(ordered_params)

    ordered_params.join
  end
end