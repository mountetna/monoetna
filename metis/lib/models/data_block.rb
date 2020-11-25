require 'fileutils'

class Metis
  class DataBlock < Sequel::Model
    plugin :timestamps, update_on_create: true

    one_to_many :files

    TEMP_PREFIX='temp-'
    TEMP_MATCH=/^#{TEMP_PREFIX}/

    def self.create_from(file_name, location, copy=false)
      # we don't know the true md5 so we use a random value
      data_block = create(
        md5_hash: "#{TEMP_PREFIX}#{Metis.instance.sign.uid}",
        description: "Originally for #{file_name}",
      )

      data_block.set_file_data(location, copy)

      return data_block
    end

    def set_file_data(file_path, copy=false)
      # Rename the existing file.
      if copy
        ::FileUtils.copy(
          file_path,
          location
        )
      else
        ::File.rename(
          file_path,
          location
        )
      end
    end

    def temp_hash?
      md5_hash =~ TEMP_MATCH
    end

    def compute_hash!
      if has_data? && temp_hash?
        md5_hash = Metis::File.md5(location)

        existing_block = Metis::DataBlock.where(md5_hash: md5_hash).first

        if existing_block
          # Point the files to the old block
          Metis::File.where(
            data_block_id: id
          ).update(
            data_block_id: existing_block.id
          )

          # destroy the redundant file
          ::File.delete(location)

          # destroy this redundant record
          destroy
          return
        end

        old_location = location

        update(md5_hash: md5_hash)

        new_location = location

        # Actually move the file
        ::File.rename(
          old_location,
          new_location
        )
      end
    end

    def backup!
      return if temp_hash? || archive_id

      Metis.instance.archiver.archive(self)
    end

    def actual_size
      has_data? ? ::File.size(location) : nil
    end

    def has_data?
      ::File.exists?(location)
    end

    def location
      directory = ::File.join(
        Metis.instance.config(:data_path),
        'data_blocks',
        md5_hash[0],
        md5_hash[1]
      )

      FileUtils.mkdir_p(directory) unless ::File.directory?(directory)

      ::File.expand_path(
        ::File.join(
          directory,
          md5_hash
        )
      )
    end

    def remove!
      if !removed
        delete_block!
        update(removed: true, updated_at: DateTime.now)
        Metis.instance.archiver.delete(self) if archive_id
      end
    end

    private

    def delete_block!
      if ::File.exists?(location)
        ::File.delete(location)
      end
    end
  end
end
