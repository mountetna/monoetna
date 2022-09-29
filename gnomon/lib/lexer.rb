class Gnomon
  class Lexer < RLTK::Lexer

    EXAMPLE_GRAMMAR_JSON = {
      tokens: {
          PROJ: { 'XMPL': 'The Example Project' },
          PROJECT: { 'The Example Project': 'The Example Project' },
          COHORT: { 'Patient': 'Patient group', 'Control': 'Control group' },
          COH: { 'P': 'Patient group', 'C': 'Control group' }
          TIM: { 'T': '# Timepoint' },
          DNA: { 'DNA' : 'Exome sequencing' },
          RNA: { 'RNA': 'Bulk RNASeq' },
          BSP: { 'WBC': 'Whole Blood Cells', 'FFP': 'FFPE slice' },
          SEP: { '-': '# Separator' }
      },
      aliases: {
          PROJ: 'PROJECT',
          COH: 'COHORT'
      },
      rules: {
          project: 'PROJECT',
          cohort: 'COHORT',
          subject: 'PROJ SEP COH .n',
          timepoint: '.subject SEP TIM .n',
          rna_seq: '.timepoint SEP BSP SEP RNA .n',
          dna_seq: '.timepoint SEP BSP SEP DNA .n'
      }
  }

    def initialize(grammar)
      @grammar = grammar
      rules_from_grammar
      super()
    end

    private

    def rules_from_grammar
      return if @grammar.nil?

      grammar.rules.each do |rule|
        # Test if inherits?
        self.class.rule()
      end
    end
  end
end