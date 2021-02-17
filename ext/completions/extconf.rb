# Rubygem doesn't provide a great hook for executing completion or shared data on installing a gem by default,
# but it does provide an extension installation hook in the form of extconf.rb.
# We are 'hacking' this behavior to do a user installation of completion files on the system
require 'mkmf'

# Create a dummy extension file. Without this RubyGems would abort the
# installation process
FileUtils.touch(File.join(Dir.pwd, 'wat.' + RbConfig::CONFIG['DLEXT']))

# Drop the completion file in the user's home directory so they can wire this based on their system.
Dir.chdir(File.expand_path('..', __FILE__)) do
  puts `cp ../../etna.completion ~/etna.completion`
  puts "etna.completion has been copied to your home directory.  Source it from your bashrc file to add completions for the etna command."
end

# This is normally set by calling create_makefile() but we don't need that
# method since we'll provide a dummy Makefile. Without setting this value
# RubyGems will abort the installation.
$makefile_created = true
create_makefile('ext')
