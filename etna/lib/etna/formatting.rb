module Etna
  class Formatting
    SIZE_UNITS = %W(B KiB MiB GiB TiB).freeze

    def self.as_size(number)
      return "unknown" if number.nil?

      if number.to_i < 1024
        exponent = 0
      else
        exponent = [(Math.log(number) / Math.log(1024)).to_i, SIZE_UNITS.size].min
        number  /= 1024 ** exponent
      end

      "#{number} #{SIZE_UNITS[exponent]}"
    end
  end
end