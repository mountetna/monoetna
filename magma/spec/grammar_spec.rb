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
          COH: {
            label: 'Cohort',
            values: {
              HS: 'Homo Sapien',
              MS: 'More Sapien'
            }
          },
          TMP: {
            label: 'Timepoint',
            values: {
              D: 'Day',
              DN: 'Day Negative'
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
        first: 'TOK SEP .n'
      }

      grammar = create(:grammar, project_name: 'labors', version_number: 1, config: @config)

      expect(grammar.tokens.length).to eq(4)
      expect(grammar.rules.length).to eq(@config[:rules].length)

      expect(grammar.rules['first'].regex).to eq(/^(SCC|SCG)(-)\d+$/)
      expect { grammar.model_name("I-am-a-little-tea-pot") }.to raise_error(Magma::Grammar::UnrecognizedIdentifierError)
      expect(grammar.model_name("SCG-1")).to eq('first')
    end

    it 'with numeric incrementor' do
      @config[:rules] = {
        first: 'COH',
        second: '.first SEP TOK .n',
      }

      grammar = create(:grammar, project_name: 'labors', version_number: 1, config: @config)

      expect(grammar.tokens.length).to eq(4)
      expect(grammar.rules.length).to eq(@config[:rules].length)

      expect(grammar.rules['second'].regex).to eq(/^(HS|MS)(-)(SCC|SCG)\d+$/)
      expect(grammar.model_name("MS-SCC2")).to eq('second')
    end

    it 'composite that build upon other rules' do
      @config[:rules] = {
        first: 'COH',
        second: '.first SEP TMP .n',
        third: '.second SEP TOK .n'
      }

      grammar = create(:grammar, project_name: 'labors', version_number: 1, config: @config)

      expect(grammar.tokens.length).to eq(4)
      expect(grammar.rules.length).to eq(@config[:rules].length)

      expect(grammar.rules['third'].regex).to eq(/^(HS|MS)(-)(D|DN)\d+(-)(SCC|SCG)\d+$/)
      expect(grammar.model_name("HS-D1-SCG10")).to eq('third')
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
            },
          },
          'SEP' => {
            'label' => 'Separator',
            'values' => {
              '-' => '-'
            }
          },
          'COH' => {
            'label' => 'Cohort',
            'values' => {
              'HS' => 'Homo Sapien',
              'MS' => 'More Sapien'
            }
          },
          'TMP' => {
            'label' => 'Timepoint',
            'values' => {
              'D' => 'Day',
              'DN' => 'Day Negative'
            }
          },
        },
        'rules' => {
          'first' => 'COH',
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
            'first' => 'TOK SEP .n',
            'second' => 'TOK SEP .n'
          }
        }

        validator = Magma::Grammar::Validation.new(config)

        expect(validator.valid?).to eq(false)
        expect(validator.errors).to eq(["Duplicate rule definition exists: [\"TOK SEP .n\"]"])
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
            'first' => 'TOK SEP .n',
            'second' => '.first'
          }
        }

        validator = Magma::Grammar::Validation.new(config)

        expect(validator.valid?).to eq(false)
        expect(validator.errors).to eq(["Duplicate rule definition exists: [\"TOK SEP .n\"]"])
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

      it 'uses a token multiple times' do
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
            'first' => 'TOK SEP TOK'
          }
        }

        validator = Magma::Grammar::Validation.new(config)

        expect(validator.valid?).to eq(false)
        expect(validator.errors).to eq(["Rule \"first\" contains duplicate tokens."])
      end

      it 'uses a token multiple times across rules' do
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
            'first' => 'TOK',
            'second' => '.first SEP TOK'
          }
        }

        validator = Magma::Grammar::Validation.new(config)

        expect(validator.valid?).to eq(false)
        expect(validator.errors).to eq(["Rule \"second\" contains duplicate tokens."])
      end

      it 'uses a numeric incrementor in the middle of a rule' do
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
            'first' => 'TOK .n SEP',
          }
        }

        validator = Magma::Grammar::Validation.new(config)

        expect(validator.valid?).to eq(false)
        expect(validator.errors).to eq(["Rule \"first\" can only use the numeric increment \".n\" at the end."])
      end
    end
  end
end