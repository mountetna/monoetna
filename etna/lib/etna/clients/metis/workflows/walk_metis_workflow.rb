module Etna
  module Clients
    class Metis
      class WalkMetisWorkflow < Struct.new(:metis_client, :project_name,
          :bucket_name, :logger, :root_dir, keyword_init: true)
        def each(&block)
          q = [self.root_dir]

          while (n = q.pop)
            resp = metis_client.list_folder(Etna::Clients::Metis::ListFolderRequest.new(
              project_name: project_name,
              bucket_name: bucket_name,
              folder_path: n
            ))

            resp.files.all.sort_by { |f| f.file_path[self.root_dir.length..-1] }.each do |file|
              yield [file, file.file_path[self.root_dir.length..-1]]
            end

            resp.folders.all.sort_by { |f| f.folder_path[self.root_dir.length..-1] }.each do |f|
              yield [f, f.folder_path[self.root_dir.length..-1]]
              q << f.folder_path
            end
          end
        end
      end
    end
  end
end