describe Magma::Grammar do

  context 'builds rules' do
    before(:each) do
      @config = {
        tokens: {
          TOK: {
            label: 'Token',
            values: {
              SCC: 'Single-cell CITEseq',
              SCG: 'Single-cell GEX'
            }
          },
          SEP: {
            label: 'Separator',
            values: {
              '-': '-'
            }
          }
        }
      }
    end

    it 'simple' do
      @config[:rules] = {
        first: 'TOK SEP TOK'
      }

      grammar = create(:grammar, project_name: 'labors', version_number: 1, config: @config)

      expect(grammar.tokens.length).to eq(2)
      expect(grammar.rules.length).to eq(@config[:rules].length)

      expect(grammar.rules['first'].regex).to eq(/^(SCC|SCG)(-)(SCC|SCG)$/)
      expect { grammar.model_name("I-am-a-little-tea-pot") }.to raise_error(Magma::Grammar::UnrecognizedIdentifierError)
      expect(grammar.model_name("SCG-SCC")).to eq('first')
    end

    it 'with numeric incrementor' do
      @config[:rules] = {
        first: 'TOK SEP TOK',
        second: '.first SEP TOK .n',
      }

      grammar = create(:grammar, project_name: 'labors', version_number: 1, config: @config)

      expect(grammar.tokens.length).to eq(2)
      expect(grammar.rules.length).to eq(@config[:rules].length)

      expect(grammar.rules['second'].regex).to eq(/^(SCC|SCG)(-)(SCC|SCG)(-)(SCC|SCG)\d+$/)
      expect(grammar.model_name("SCG-SCC-SCC2")).to eq('second')
    end

    it 'composite that build upon other rules' do
      @config[:rules] = {
        first: 'TOK SEP TOK',
        second: '.first SEP TOK .n',
        third: '.second SEP TOK .n'
      }

      grammar = create(:grammar, project_name: 'labors', version_number: 1, config: @config)

      expect(grammar.tokens.length).to eq(2)
      expect(grammar.rules.length).to eq(@config[:rules].length)

      expect(grammar.rules['third'].regex).to eq(/^(SCC|SCG)(-)(SCC|SCG)(-)(SCC|SCG)\d+(-)(SCC|SCG)\d+$/)
      expect(grammar.model_name("SCG-SCC-SCC1-SCG10")).to eq('third')
    end
  end

  describe Magma::Grammar::Validation do
    it 'throws no errors on valid rules' do
      config = {
        'tokens' => {
          'TOK' => {
            'label' => 'Token',
            'values' => {
              'SCC' => 'Single-cell CITEseq',
              'SCG' => 'Single-cell GEX'
            }
          },
          'SEP' => {
            'label' => 'Separator',
            'values' => {
              '-' => '-'
            }
          }
        },
        'rules' => {
          'first' => 'TOK SEP TOK',
          'second' => '.first SEP TOK .n'
        }
      }

      validator = Magma::Grammar::Validation.new(config)

      expect(validator.valid?).to eq(true)
      expect(validator.errors).to eq([])
    end

    context 'complains if rule' do
      it 'uses unspecified token' do
        config = {
          'tokens' => {
            'TOK' => {
              'label' => 'Token',
              'values' => {
                'SCC' => 'Single-cell CITEseq',
                'SCG' => 'Single-cell GEX'
              }
            },
            'SEP' => {
              'label' => 'Separator',
              'values' => {
                '-' => '-'
              }
            }
          },
          'rules' => {
            'first' => 'ABC SEP TOK'
          }
        }

        validator = Magma::Grammar::Validation.new(config)

        expect(validator.valid?).to eq(false)
        expect(validator.errors).to eq(["Unknown tokens used: [\"ABC\"]"])
      end

      it 'uses unspecified other rule' do
        config = {
          'tokens' => {
            'TOK' => {
              'label' => 'Token',
              'values' => {
                'SCC' => 'Single-cell CITEseq',
                'SCG' => 'Single-cell GEX'
              }
            },
            'SEP' => {
              'label' => 'Separator',
              'values' => {
                '-' => '-'
              }
            }
          },
          'rules' => {
            'first' => '.imaginary SEP TOK'
          }
        }

        validator = Magma::Grammar::Validation.new(config)

        expect(validator.valid?).to eq(false)
        expect(validator.errors).to eq(["Unknown tokens used: [\".imaginary\"]"])
      end

      it 'is duplicated' do
        config = {
          'tokens' => {
            'TOK' => {
              'label' => 'Token',
              'values' => {
                'SCC' => 'Single-cell CITEseq',
                'SCG' => 'Single-cell GEX'
              }
            },
            'SEP' => {
              'label' => 'Separator',
              'values' => {
                '-' => '-'
              }
            }
          },
          'rules' => {
            'first' => 'TOK SEP TOK',
            'second' => 'TOK SEP TOK'
          }
        }

        validator = Magma::Grammar::Validation.new(config)

        expect(validator.valid?).to eq(false)
        expect(validator.errors).to eq(["Duplicate rule definition exists: [\"TOK SEP TOK\"]"])
      end

      it 'is duplicated using a rule name' do
        config = {
          'tokens' => {
            'TOK' => {
              'label' => 'Token',
              'values' => {
                'SCC' => 'Single-cell CITEseq',
                'SCG' => 'Single-cell GEX'
              }
            },
            'SEP' => {
              'label' => 'Separator',
              'values' => {
                '-' => '-'
              }
            }
          },
          'rules' => {
            'first' => 'TOK SEP TOK',
            'second' => '.first'
          }
        }

        validator = Magma::Grammar::Validation.new(config)

        expect(validator.valid?).to eq(false)
        expect(validator.errors).to eq(["Duplicate rule definition exists: [\"TOK SEP TOK\"]"])
      end

      it 'is recursive' do
        config = {
          'tokens' => {
            'TOK' => {
              'label' => 'Token',
              'values' => {
                'SCC' => 'Single-cell CITEseq',
                'SCG' => 'Single-cell GEX'
              }
            },
            'SEP' => {
              'label' => 'Separator',
              'values' => {
                '-' => '-'
              }
            }
          },
          'rules' => {
            'first' => '.second SEP TOK',
            'second' => '.first'
          }
        }

        validator = Magma::Grammar::Validation.new(config)

        expect(validator.valid?).to eq(false)
        expect(validator.errors).to eq(["Rule \"second\" may be recursive! It's token \".first\" appears to lead to circular logic.", "Rule \"first\" may be recursive! It's token \".second\" appears to lead to circular logic."])
      end

      it 'is blank' do
        config = {
          'tokens' => {
            'TOK' => {
              'label' => 'Token',
              'values' => {
                'SCC' => 'Single-cell CITEseq',
                'SCG' => 'Single-cell GEX'
              }
            },
            'SEP' => {
              'label' => 'Separator',
              'values' => {
                '-' => '-'
              }
            }
          },
          'rules' => {
            'first' => ''
          }
        }

        validator = Magma::Grammar::Validation.new(config)

        expect(validator.valid?).to eq(false)
        expect(validator.errors).to eq(["Rules cannot contain only whitespace"])
      end

      it 'is self-referential' do
        config = {
          'tokens' => {
            'TOK' => {
              'label' => 'Token',
              'values' => {
                'SCC' => 'Single-cell CITEseq',
                'SCG' => 'Single-cell GEX'
              }
            },
            'SEP' => {
              'label' => 'Separator',
              'values' => {
                '-' => '-'
              }
            }
          },
          'rules' => {
            'first' => '.first SEP TOK'
          }
        }

        validator = Magma::Grammar::Validation.new(config)

        expect(validator.valid?).to eq(false)
        expect(validator.errors).to eq(["Rule \"first\" may be recursive! It's token \".first\" appears to lead to circular logic."])
      end
    end
  end
end