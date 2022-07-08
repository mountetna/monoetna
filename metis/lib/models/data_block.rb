require "fileutils"
require "tempfile"
require 'securerandom'

class Metis

  class ThumbnailError < StandardError
  end

  class ThumbnailNotExistError < ThumbnailError
  end

  class DataBlock < Sequel::Model
    plugin :timestamps, update_on_create: true

    one_to_many :files

    TEMP_PREFIX = "temp-"
    TEMP_MATCH = /^#{TEMP_PREFIX}/

    def self.create_from(file_name, location, copy: false)
      # we don't know the true md5 so we use a random value
      size = ::File.size(location)
      data_block = create(
        md5_hash: "#{TEMP_PREFIX}#{Metis.instance.sign.uid}",
        description: "Originally for #{file_name}",
        size: size,
      )

      data_block.set_file_data(location, copy)

      return data_block
    end

    def set_file_data(file_path, copy = false)
      if copy
        ::FileUtils.copy(
          file_path,
          location
        )
      else
        # We're running into an issue with Metis uploads where the upload
        #   cannot be renamed when creating the data_block.
        #   > Errno::EBUSY: Device or resource busy @ rb_file_s_rename -
        #   > File "/var/www/metis/lib/models/data_block.rb", line 32, in rename
        #   > File "/var/www/metis/lib/models/data_block.rb", line 19, in create_from
        #   > File "/var/www/metis/lib/models/upload.rb", line 55, in finish!
        #   > File "/var/www/metis/lib/server/controllers/upload_controller.rb", line 150, in complete_upload
        #   > File "/var/www/metis/lib/server/controllers/upload_controller.rb", line 117, in upload_blob
        #   > File "/var/www/metis/lib/models/data_block.rb", line 32, in set_file_data
        # I suspect it's a BeeGFS issue, specifically with how
        #   files get closed. Note the instructions for the config flag `tuneEarlyCloseResponse`:
        #   > # Request close responses from the metadata server before the file is fully closed.
        #   > # This may improve close() performance, but closed files may be accounted as
        #   > # open for a short time after close() has returned. Files accounted as open
        #   > # cannot be moved.
        #   > # Default: false
        # So it seems like files may be kept open for some amount
        #   of time after the final blob is appended, which will
        #   prevent it from being moved.
        # Adding a delay (super hacky?!?) will hopefully give the
        #   file system time to close the file. Hopefully this does not
        #   significantly slow down uploads.
        sleep(0.01)
        ::File.rename(
          file_path,
          location
        )
      end
    end

    def temp_hash?
      md5_hash =~ TEMP_MATCH
    end

    def compute_hash!(&test_hook)
      if temp_hash? && has_data?
        md5_hash = Metis::File.md5(location)
        existing_block = Metis::DataBlock.where(md5_hash: md5_hash).first

        if existing_block
          # Point the files to the old block
          Metis::File.where(
            data_block_id: id,
          ).update(
            data_block_id: existing_block.id,
          )

          yield [:temp_delete] if block_given?
          # destroy the redundant file
          ::File.delete(location) if has_data?

          # destroy this redundant record
          yield [:temp_destroy] if block_given?
          destroy
          return
        else
          # Create a stub in position without removing the exist file, to ensure that
          # if the database transaction fails, we still have a way of finding the file backing
          # this datablock.
          yield [:sym, location, file_location(md5_hash)] if block_given?
          begin
            ::File.symlink(location, file_location(md5_hash))
          rescue ::Errno::EEXIST => e
            # Some other process may have already prepared the file by the same hash,
          end

          yield [:new_update] if block_given?
          update(md5_hash: md5_hash)
        end
      end

      if ::File.symlink?(location) && !temp_hash?
        original_file = ::File.readlink(location)
        ::File.rename(original_file, location)
      end

      unless has_data?
        yield [:update_removed] if block_given?
        update(removed: true)
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
      file_location(md5_hash)
    end

    def named_linked_location(file)
      ::File.join(
        Metis.instance.config(:data_path),
        "named",
        file.bucket.name,
        file.file_path
      )
    end

    def file_location(hash)
      directory = ::File.join(
        Metis.instance.config(:data_path),
        "data_blocks",
        hash[0],
        hash[1]
      )

      FileUtils.mkdir_p(directory) unless ::File.directory?(directory)

      ::File.expand_path(
        ::File.join(
          directory,
          hash
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

    def thumbnail
      raise ThumbnailNotExistError.new("Thumbnail does not exist for data_block #{md5_hash}") unless thumbnail_in_cache?

      cached_thumbnail
    end

    def thumbnail_in_cache?
      ::File.exist?(thumbnail_location)
    end

    def generate_thumbnail
      begin
        update(has_thumbnail: false)
        return
      end unless has_data?

      begin
        source = Vips::SourceCustom.new
        ::File.open(location, "rb") do |handle|
          source.on_read { |length| handle.read(length) }
          source.on_seek { |offset, whence| handle.seek(offset, whence) }

          th = Vips::Image.thumbnail_source(source, 640, height: 480)
          th.write_to_file(thumbnail_location)
        end
        update(has_thumbnail: true)
      rescue Vips::Error => e
        update(has_thumbnail: false)
        Metis.instance.logger.error(e.message)
      ensure
        refresh
      end
    end

    private

    def delete_block!
      if ::File.exists?(location)
        ::File.delete(location)
      end
    end

    def cached_thumbnail
      ::File.read(thumbnail_location)
    end

    def thumbnail_location
      directory = ::File.join(
        Metis.instance.config(:thumbnail_path),
        md5_hash[0],
        md5_hash[1],
      )

      FileUtils.mkdir_p(directory) unless ::File.directory?(directory)

      ::File.expand_path(
        ::File.join(
          directory,
          "#{::File.basename(location)}.png"
        )
      )
    end
  end
end
