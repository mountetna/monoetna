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
      update(md5_hash: Metis::File.md5(location)) if has_data? && temp_hash?
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
      ::File.expand_path(
        ::File.join(
          Metis.instance.config(:data_path),
          'data_blocks',
          md5_hash
        )
      )
    end
  end
end
