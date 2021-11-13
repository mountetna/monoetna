module Etna
  module Clients
    class Metis
      class WalkMetisDiffWorkflow < Struct.new(:left_walker, :right_walker, keyword_init: true)
        def each(&block)
          left_enum = self.left_walker.to_enum
          right_enum = self.right_walker.to_enum

          l, l_path = next_or_nil(left_enum)
          r, r_path = next_or_nil(right_enum)

          while l && r
            if l_path == r_path
              yield [compare_file_or_folders(l, r), l, r]

              l, l_path = next_or_nil(left_enum)
              r, r_path = next_or_nil(right_enum)
            elsif l_path < r_path
              yield [:left_unique, l, nil]
              l, l_path = next_or_nil(left_enum)
            else
              yield [:right_unique, nil, r]
              r, r_path = next_or_nil(right_enum)
            end
          end

          while l
            yield [:left_unique, l, nil]
            l, l_path = next_or_nil(left_enum)
          end

          while r
            yield [:right_unique, nil, r]
            r, r_path = next_or_nil(right_enum)
          end
        end

        def next_or_nil(enum)
          enum.next
        rescue StopIteration
          [nil, nil]
        end

        def compare_file_or_folders(l, r)
          if l.is_a?(Etna::Clients::Metis::Folder)
            if r.is_a?(Etna::Clients::Metis::Folder)
              return compare_file_or_folder_age(l, r)
            end

            return :left_is_folder
          end

          if r.is_a?(Etna::Clients::Metis::Folder)
            return :right_is_folder
          end


          if l.file_hash.nil? || r.file_hash.nil?
            return :unknown
          end

          if l.file_hash == r.file_hash
            return :equal
          end

          return compare_file_or_folder_age(l, r)
        end

        def compare_file_or_folder_age(l, r)
          if l.updated_at.nil?
            return :unknown
          end

          if r.updated_at.nil?
            return :unknown
          end

          if l.updated_at < r.updated_at
            return :left_older
          elsif l.updated_at > r.updated_at
            return :right_older
          else
            return :equal
          end
        end
      end
    end
  end
end
