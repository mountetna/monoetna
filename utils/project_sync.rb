#!/home/developer/.rbenv/shims/ruby

# project_sync.rb
# This script will read the projects out of Janus and make sure that there are
# corresponding folders to house the files.

require 'net/http'
require 'json'
require 'fileutils'
require 'openssl'
require '../server/conf'
require '../server/secrets'

module ProjectSync

  def ProjectSync.get_groups()

    url = Conf::JANUS_ADDR
    url = url + '/get-groups'
    data = {

      :token=> ARGV[0], 
      :app_key=> Secrets::APP_KEY
    }

    groups = make_request(url, data)

    return groups
  end

  def ProjectSync.check_groups(groups)

    groups = JSON.parse(groups)

    if !groups.key?('success') || !groups.key?('groups')

      puts ProjectSync.send_server_error()
      return 0
    end

    if !groups['success']

      puts ProjectSync.send_server_error()
      return 0
    end

    return groups['groups']
  end
  
  def ProjectSync.get_projects()

    url = Conf::JANUS_ADDR
    url = url + '/get-projects'
    data = {

      :token=> ARGV[0], 
      :app_key=> Secrets::APP_KEY 
    }
    projects = make_request(url, data)
    return projects
  end

  def ProjectSync.check_projects(projects)

    projects = JSON.parse(projects)
    if !projects.key?('success') || !projects.key?('projects')

      puts ProjectSync.send_server_error()
      return 0
    end 

    if !projects['success']

      puts ProjectSync.send_bad_request()
      return 0
    end

    return projects['projects']
  end

  def ProjectSync.set_group_folders(structure)

    # Create any directories on the file system that are new.
    structure.each do |group_key, group_projects|

      dir_name = File.join('/data', group_key.to_s())
      if !File.directory?(dir_name)

        FileUtils.mkdir_p(dir_name)
      end

      ProjectSync.set_project_folders(group_key, group_projects)
    end
  end

  def ProjectSync.set_project_folders(group_key, group_projects)

    group_projects.each do |project_key, file_names|

      sub_dir_name = File.join('/data', group_key.to_s(), project_key.to_s())
      if !File.directory?(sub_dir_name)

        FileUtils.mkdir_p(sub_dir_name)
      end
    end
  end

  def ProjectSync.normalize_name(nm)

    nm.gsub!('/ ', '_')
    nm.gsub!('(', '')
    nm.gsub!(')', '')
    nm.gsub!(' ', '_')
    nm.gsub!('/', '_')
    nm.gsub!('-', '_')
    return nm
  end

  def ProjectSync.extract_structure(groups, projects)

    structure = {}
    groups.each do |group|

      grp_nm = ProjectSync.normalize_name(group['group_name'])
      structure[grp_nm] = {}
    end

    # Create the directory structure from the database.
    projects.each do |project|

      grp_nm = ProjectSync.normalize_name(project['group_name'])
      if !structure.key?(grp_nm)

        # alert that a main group is missing.
        # it is present in the projects but no in the group table.
      else

        prjk_nm = ProjectSync.normalize_name(project['project_name'])
        structure[grp_nm][prjk_nm] = []
      end
    end

    return structure
  end

  def ProjectSync.extract_new_structure()

    new_structure = {}
    Dir.chdir('/data')
    Dir.glob('*').select {|dir_file| 

      if File.directory?(dir_file)

        # puts dir_file
        new_structure[dir_file] = {}
      end
    }

    return new_structure
  end

  def ProjectSync.make_request(url, data)

    begin

      uri = URI.parse(url)
      https_conn = Net::HTTP.new(uri.host, uri.port)
      https_conn.use_ssl = true
      https_conn.verify_mode = OpenSSL::SSL::VERIFY_PEER

      request = Net::HTTP::Post.new(uri.path)
      request.set_form_data(data)

      response = https_conn.request(request)
      return response.body
    rescue Timeout::Error, 
           Errno::EINVAL, 
           Errno::ECONNRESET, 
           EOFError, 
           Net::HTTPBadResponse, 
           Net::HTTPHeaderSyntaxError, 
           Net::ProtocolError => error

      return ProjectSync::send_server_error()
    end
  end

  def ProjectSync.send_server_error()

    return 'There was a server error.'
  end

  def ProjectSync.print_help

    puts ''
    puts ' desc: Generates a tree that represents the metis file/directory structure.'
    puts 'usage: ./project_sync.rb [janus token]'
    puts ''
  end
end

# If the argument is 'nil', or if there are to many arguments, or if the
# argument is an empty string, then bail.
if ARGV[0].nil? || ARGV.length != 1 || ARGV[0].to_s.length == 0

  ProjectSync::print_help()
  exit 1
end

# Get and check the groups from Janus
groups = ProjectSync::get_groups()
groups = ProjectSync::check_groups(groups)

if groups == 0

  exit 1
end

# Get and check the projects from Janus.
projects = ProjectSync::get_projects()
projects = ProjectSync::check_projects(projects)

if projects == 0

  exit 1
end

# Extract the structure of the directory from the database.
structure = ProjectSync::extract_structure(groups, projects)
# Create any directories on the file system that are new.

ProjectSync::set_group_folders(structure)