require_relative 'revision'

class Metis
    class CopyRevision < Revision
        def self.create_from_parts(params)
            Metis::CopyRevision.new({
              source: Metis::Path.path_from_parts(
                params[:source][:project_name],
                params[:source][:bucket_name],
                params[:source][:file_path]),
              dest: Metis::Path.path_from_parts(
                params[:dest][:project_name],
                params[:dest][:bucket_name],
                params[:dest][:file_path])
            })
        end

        def validate (user_authorized_bucket_names)
          super(user_authorized_bucket_names)

          @errors.push("Invalid dest path: #{@dest.path}") unless @dest.valid?

          if @dest.valid?
            @errors.push("Invalid dest bucket: #{@dest.bucket_name}") unless @dest.bucket

            if @dest.bucket
              @errors.push(
                "Forbidden: no access to dest bucket #{@dest.bucket.name}"
              ) unless user_authorized_bucket_names.include? @dest.bucket.name

              dest_folder = @dest.folder

              if @dest.folder_path
                @errors.push("Invalid dest folder: #{@dest.folder_path}") unless dest_folder
                @errors.push("Dest folder #{dest_folder.folder_name} is read-only") if dest_folder&.read_only?
              end

              dest_file = @dest.file

              if dest_file
                @errors.push("Dest file #{@dest.file_name} is read-only") if dest_file.read_only?
              end

              @errors.push(
                "Cannot copy over existing folder #{@dest.path}"
              ) if  Metis::Folder.exists?(@dest.file_name, @dest.bucket, dest_folder)

            end
          end
        end

        def bucket_names
          # This is kind of weird, but we need the ability to grab
          #   all relevant bucket names, even before validation of
          #   the CopyRevision (in file_controller).
          # So, if @dest doesn't exist, we only return the source
          #   bucket name.
          if @dest.valid?
            return [@source.bucket_name, @dest.bucket_name]
          end
          return [@source.bucket_name]
        end
    end
end