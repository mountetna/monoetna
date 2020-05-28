require_relative 'path'

class Metis
    class PathWithObjects

      # This class is for convenience to store the query-able
      #     objects, so in the bulk operations we
      #     can minimize database calls.
      attr_reader :mpath
      attr_accessor :bucket, :folder, :file
      def initialize(mpath)
        @mpath = Metis::Path.new(mpath)
      end
    end
  end
