require 'net/http/post/multipart'
require 'singleton'
require_relative '../base_client'
require_relative '../../client'
require_relative './models'

module Etna
  module Clients
    class Gnomon < Etna::Clients::BaseClient
      def project_rules(project_name)
        json = nil
        @etna_client.get("/gnomon/#{project_name}/rules") do |res|
          json = JSON.parse(res.body)
        end

        ProjectRulesResponse.new(json)
      end
    end
  end
end
