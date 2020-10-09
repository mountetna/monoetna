module RegexDsl
  def regex_chain(*regs)
    flags = 0
    src = regs.map do |reg|
      case reg
      when Regexp
        [reg.source]
      when String
        [reg]
      when Integer
        flags |= reg
        []
      else
        raise "Got value #{reg}, expected regex"
      end
    end.inject([], &:+).join('\.')

    Regexp.new(src, flags)
  end

  def regex_named_match(name, reg)
    source =
      case reg
      when Regexp
        reg.source
      when String
        reg
      else
        raise "Got value #{reg}, expected regex"
      end
    Regexp.new("(?<#{name}>#{source})")
  end

  def regex_options(options)
    /(?:#{options.join('|')})/
  end
end