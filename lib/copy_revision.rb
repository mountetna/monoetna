require_relative 'revision'

class Metis
  class CopyRevision < Revision
    def validate
      @errors = []
      validate_source(@source)
      validate_dest(@dest)
    end

    def bucket_names
      # This is kind of weird, but we need the ability to grab
      #   all relevant bucket names, even before validation of
      #   the CopyRevision.
      # So, if @dest doesn't exist, we don't return it
      source_bucket_names = super
      if @dest.mpath.valid?
        return source_bucket_names.push(@dest.mpath.bucket_name)
      end
      return source_bucket_names
    end

    def mpaths
      source_mpaths = super
      source_mpaths.push(@dest.mpath) if @dest.mpath.valid?
      return source_mpaths
    end

    def revise!
      raise 'Invalid revision, cannot revise!' unless valid?
      raise 'Cannot revise without a user' unless @user

      Metis::File.copy({
        project_name: @source.mpath.project_name,
        source_file: @source.file,
        dest_file_path: @dest.mpath.file_path,
        dest_bucket_name: @dest.mpath.bucket_name,
        user: @user
      })
    end

    def to_hash
      {
        dest: @dest.mpath.path,
        source: @source.mpath.path,
        errors: @errors
      }
    end
  end
end