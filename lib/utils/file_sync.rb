#!/home/developer/.rbenv/shims/ruby

# This script will read a directory and the Metis database and try to make a
# match. It will match names and md5 hashes.

# 1 get a list of all the files in a project directory. Ignore subdirectories

require('digest/md5')
require('json')
require('sequel')

require('../server/secrets')
require('../server/service/postgres_service.rb')
PostgresService::connect
require('../server/models/file_model.rb')

module FileSync
  def FileSync.get_project_files
    if File.directory?("/data1/#{ARGV[0]}")
      Dir.entries("/data1/#{ARGV[0]}").select do |f|
        next(false) if File.directory?(f)       # Ignore if it's a directory.
        next(false) if f == 'hashes.json'       # Ignore if it's a hash summary.
        next(false) if File.extname(f)=='.part' # Ignore if it's a partial file.
        next(true)
      end
    end
  end

  def FileSync.hash_file(file_name, size)
    md5 = Digest::MD5.new
    File.open("/data1/#{ARGV[0]}/#{file_name}").each(nil, size) do |chunk|
      md5.update(chunk)
    end
    md5.hexdigest
  end

  def FileSync.save_hashes(hashes)
    File.open("/data1/#{ARGV[0]}/hashes.json", 'w') do |f|
      f.write(JSON.pretty_generate(hashes))
    end
  end

  def FileSync.get_database_entries
    PostgresService::get_files_by_project_name(ARGV[0])
  end

  # Also need to check for a hash mismatch
  def FileSync.check_for_missing_files(db_entries, file_hashes)
    db_entries.select { |entry| !file_hashes.key?(entry[:file_name]) }
  end

  def FileSync.check_for_extras_files(file_hashes, db_entries)
    file_hashes.select do |file_name, file_hash|
      entries = db_entries.select do |entry|
        file_name == entry[:file_name] ? true : false
      end
      entries.length > 0 ? false : true
    end
  end

  def FileSync.add_extra_files(extra_files, group_name)
    timestamp = Time::now
    extra_files.each do |key, value|
      params = {
        'original_name'=> key,
        'file_name'=> key,
        'file_size'=> File.size("/data1/#{ARGV[0]}/#{key}"),
        'group_name'=> group_name,
        'project_name'=> ARGV[0],
        'project_name_full'=>'immunoprofiler initiative',
        'start_upload'=> timestamp,
        'finish_upload'=>  timestamp,
        'user_email'=> 'jason.cater@ucsf.edu',
        'hashing_algorithm'=> 'MD5',
        'hash'=> value
      }
      PostgresService::create_new_file!(params)
    end
  end

  def FileSync.prompt(*args)
    puts(*args)
    STDIN.gets.chomp
  end

  def FileSync.print_help
    puts ''
    puts ' desc: Syncs the files in a project directory with the entries in '\
'the Metis databse.'
    puts 'usage: ./file_sync.rb [project_name]'
    puts ''
  end
end

# If the argument is 'nil', or if there are to many arguments, or if the
# argument is an empty string, then bail.
if ARGV[0].nil? || ARGV.length != 1 || ARGV[0].to_s.length == 0
  FileSync::print_help
  exit 1
end

# Get all the file names from the project folder.
file_names = FileSync::get_project_files

# Hash the files from the project and save to an object. We are also setting the
# blob/chunk size to hash at a time to 1MiB.
file_hashes = {}
file_names.each do |file_name|
  file_hashes[file_name] = FileSync.hash_file(file_name, 2**20)
end

# Save the hashes to a file for later use.
FileSync.save_hashes(file_hashes)

# Extract the entries for a project 
db_entries = FileSync::get_database_entries

# Check if there are missing files. i.e. we have a db record but no file exists.
missing_files = FileSync::check_for_missing_files(db_entries, file_hashes)
puts('')
puts('missing files:')
puts('--------------------------------------------------------------------------------')
if missing_files.length == 0
  puts('[]')
else
  puts(JSON.pretty_generate(missing_files))
end
puts('')

# Check for extra files that are present in the project directory but do not
# have an entry in the db.
extra_files = FileSync::check_for_extras_files(file_hashes, db_entries)
puts('extra files:')
puts('--------------------------------------------------------------------------------')
puts(JSON.pretty_generate(extra_files))
puts('')

# Ask if we want to add the extra files.
answer = FileSync::prompt('Do you want to add the extra files to the db? (y/n)')
if answer == 'y'
  group_name = FileSync::prompt('What is the group association?')
  exit 1 if group_name == ''
  FileSync::add_extra_files(extra_files, group_name)
end
