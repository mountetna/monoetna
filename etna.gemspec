Gem::Specification.new do |spec|
  spec.name              = 'etna'
  spec.version           = '0.1.2'
  spec.summary           = 'Base classes for Mount Etna applications'
  spec.description       = 'See summary'
  spec.email             = 'Saurabh.Asthana@ucsf.edu'
  spec.homepage          = 'http://github.com/mountetna/etna'
  spec.author            = 'Saurabh Asthana'
  spec.files             = Dir['lib/**/*.rb']
  spec.platform          = Gem::Platform::RUBY # This is the default
  spec.require_paths     = [ 'lib' ]
  spec.extensions        = Dir['ext/**/extconf.rb']
  spec.add_dependency 'extlib'
  spec.add_dependency 'rack'
  spec.add_dependency 'jwt'
end
