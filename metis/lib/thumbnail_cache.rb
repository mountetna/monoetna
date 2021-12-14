require "mimemagic"
require "vips"

class Metis
  class ThumbnailNotExistError < StandardError
  end

  class ThumbnailCache
    def mimetype(file)
      MimeMagic.by_path(file.file_name)
    end

    def thumbnail(file)
      raise ThumbnailNotExistError.new("Thumbnail does not exist for file #{file.name}") unless thumbnail_in_cache?(file)

      cached_thumbnail(file)
    end

    def thumbnail_in_cache?(file)
      return false if file.data_block.nil?

      ::File.exist?(thumbnail_path(file))
    end

    def generate_thumbnail(file)
      thumbnail = Vips::Image.thumbnail(file.data_block.location, 240)
      mimetype = mimetype(file)

      case mimetype.to_s
      when /jpe?g$/
        buffer = thumbnail.jpegsave_buffer
      when "image/tiff"
        buffer = thumbnail.tiffsave_buffer
      when "image/png"
        buffer = thumbnail.pngsave_buffer
      else
        raise Etna::Error.new("Thumbnails not supported for mimetype #{mimetype}", 422)
      end

      cache_thumbnail(file, buffer)

      buffer
    end

    private

    def cache_thumbnail(original_file, thumbnail_data)
      ::File.open(thumbnail_path(original_file), "w") do |f|
        f.write(thumbnail_data)
      end
    end

    def cached_thumbnail(file)
      ::File.read(thumbnail_path(file))
    end

    def thumbnail_path(file)
      ::File.join(
        ::File.dirname(file.data_block.location),
        "th_#{::File.basename(file.data_block.location)}"
      )
    end
  end
end
