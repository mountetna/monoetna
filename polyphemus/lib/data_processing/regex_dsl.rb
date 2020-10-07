module RegexDsl
  def regex_chain(*regs)
    Regexp.new(regs.map do |reg|
      case reg
      when Regexp
        reg.source
      when String
        reg
      else
        raise "Got value #{reg}, expected regex"
      end
    end.join('\.'))
  end

  def regex_options_match(options)
    /(?:#{options.join('|')})/
  end
end