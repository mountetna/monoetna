module Enum
  def self.included(cls)
    cls.instance_eval do
      def entries
        self.constants.map { |c| self.const_get(c) }.select { |v| v.instance_of?(self) }
      end
    end
  end
end