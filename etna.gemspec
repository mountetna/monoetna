Gem::Specification.new do |spec|
  spec.name              = 'etna'
  spec.version           = '0.1.39'
  spec.summary           = 'Base classes for Mount Etna applications'
  spec.description       = 'See summary'
  spec.email             = 'Saurabh.Asthana@ucsf.edu'
  spec.homepage          = 'http://github.com/mountetna/etna'
  spec.author            = 'Saurabh Asthana'
  spec.files             = Dir['lib/**/*.rb', 'lib/etna/templates/*.json', '*.completion']
  spec.platform          = Gem::Platform::RUBY # This is the default
  spec.require_paths     = [ 'lib' ]
  spec.extensions        = Dir['ext/**/extconf.rb']
  spec.executables       << 'etna'
  spec.add_dependency 'rack'
  spec.add_dependency 'jwt'
  spec.add_dependency 'nokogiri'
  spec.add_dependency 'multipart-post'
  spec.add_dependency 'rollbar'
  spec.add_dependency 'concurrent-ruby'
  spec.add_dependency 'net-sftp', '>=3.0.0'
  spec.post_install_message = "\n\n\033[1;31mAdd 'source ~/etna.completion' to your .bashrc to get etna command line completions!\033[0m\n\n"
end
