#!/home/developer/.rbenv/shims/ruby

# This script will read the projects out of Janus and make sure that there are
# corresponding folders to house the files.

require 'net/http'
require 'json'
require 'fileutils'
require 'openssl'
require '../server/conf'
require '../server/secrets'

module ProjectSync
  def ProjectSync.get_projects
    url = "#{Conf::JANUS_ADDR}/get-projects"
    projects = make_request(url, {token: ARGV[0], app_key: Secrets::APP_KEY})
    return projects
  end

  def ProjectSync.check_projects(projects)
    projects = JSON.parse(projects)

    if !projects.key?('success') || !projects.key?('projects')
      puts ProjectSync.send_server_error
      return 0
    end 

    if !projects['success']
      puts ProjectSync.send_server_error
      return 0
    end

    return projects['projects']
  end

  def ProjectSync.set_project_folders(projects)
    projects.each do |project|
      dir_name = File.join('/data1', project['project_name'])
      FileUtils.mkdir_p(dir_name) if !File.directory?(dir_name)
    end
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
      return ProjectSync::send_server_error
    end
  end

  def ProjectSync.send_server_error
    'There was a server error.'
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
  ProjectSync::print_help
  exit 1
end

# Get and check the projects from Janus.
projects = ProjectSync::get_projects
projects = ProjectSync::check_projects(projects)
exit 1 if projects == 0
ProjectSync::set_project_folders(projects)
