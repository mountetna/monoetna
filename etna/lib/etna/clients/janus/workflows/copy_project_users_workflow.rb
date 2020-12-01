require 'base64'
require 'json'
require 'ostruct'
require_relative '../models'

module Etna
  module Clients
    class Janus
      class CopyProjectUsersWorkflow < Struct.new(:janus_client, :project_name, keyword_init: true)
        COLUMNS = [
            :name,
            :email,
            :role,
            :affiliation,
        ]

        def each_user_template_csv_row
          yield COLUMNS
          require 'oga'

          # TODO:  Replace with a real api at some point.
          parsed = Oga.parse_html(janus_client.get_project(GetProjectRequest.new(project_name: project_name)).html)

          parsed.css('.items .item').each do |item|
            name = item.at_css('.name').text.chomp
            email = item.at_css('.email').text.chomp
            role = item.at_css('.role').at_css('option[selected=""]').text.chomp
            affiliation = item.at_css('.affiliation').attribute('value').first.value

            yield row_from_columns(name: name, email: email, role: role, affiliation: affiliation)
          end
        end

        def copy_user_row(user_row)
          role = user_row[:role]
          role = 'editor' if role == 'administrator'

          janus_client.add_user(Etna::Clients::Janus::AddUserRequest.new(
              project_name: project_name,
              email: user_row[:email],
              name: user_row[:name],
              role: role,
              affiliation: user_row[:affiliation]
          ))

          if user_row[:role] == 'administrator'
            janus_client.update_permission(Etna::Clients::Janus::UpdatePermissionRequest.new(
                project_name: project_name,
                email: user_row[:email],
                role: 'administrator',
            ))
          end
        end

        def apply_user_templates_csv(io, &update_block)
          CSV.parse(io, headers: true, header_converters: :symbol).each do |row|
            copy_user_row(row)
            yield row if block_given?
          end
        end

        def write_user_templates_csv(io)
          csv = CSV.new(io)
          each_user_template_csv_row do |row|
            csv << row
          end
        end

        def row_from_columns(**columns)
          COLUMNS.map { |c| (columns[c] || '').to_s }
        end
      end
    end
  end
end

