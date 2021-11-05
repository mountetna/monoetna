import React from 'react';
import {
  render,
  fireEvent,
  waitFor,
  screen,
  within
} from '@testing-library/react';
import '@testing-library/jest-dom/extend-expect';

import {redcapSpecWrapper} from '../../spec/helpers';

import RedcapForm, {Config} from '../redcap-form';

describe('RedcapForm', () => {
  describe('renders', () => {
    it('blank for new config', async () => {
      const {asFragment} = render(
        <RedcapForm
          config={{}}
          project_name='test'
          job={{
            name: 'redcap loader',
            schema: {},
            params: {},
            secrets: {}
          }}
          update={() => {}}
        />,
        {wrapper: redcapSpecWrapper({props: {}})}
      );

      await waitFor(() => screen.getByTitle('Add model'));

      expect(asFragment()).toMatchSnapshot();

      expect(screen.queryByRole('tab')).toBeFalsy();
    });

    it('renders with existing config', async () => {
      const {asFragment} = render(
        <RedcapForm
          config={{
            testModel: {
              scripts: [
                {
                  attributes: {
                    name: 'redcap_field_name'
                  }
                }
              ],
              each: ['record'],
              invert: true,
              identifier_fields: ['non_standard_field']
            }
          }}
          project_name='test'
          job={{
            name: 'redcap loader',
            schema: {},
            params: {},
            secrets: {}
          }}
          update={() => {}}
        />,
        {
          wrapper: redcapSpecWrapper({
            props: {
              models: {
                testModel: {
                  documents: {},
                  template: {
                    attributes: {
                      name: {
                        attribute_type: 'string'
                      }
                    }
                  }
                }
              }
            }
          })
        }
      );

      await waitFor(() => screen.getByText('testModel'));

      expect(asFragment()).toMatchSnapshot();

      expect(screen.queryAllByRole('tab').length > 0).toBeTruthy();

      expect(screen.getByDisplayValue('non_standard_field')).toBeTruthy();
      expect(screen.getByDisplayValue('redcap_field_name')).toBeTruthy();
    });
  });

  describe('user can', () => {
    it('add model', async () => {
      render(
        <RedcapForm
          config={{}}
          project_name='test'
          job={{
            name: 'redcap loader',
            schema: {},
            params: {},
            secrets: {}
          }}
          update={() => {}}
        />,
        {
          wrapper: redcapSpecWrapper({
            props: {
              models: {
                testModel: {
                  documents: {},
                  template: {
                    attributes: {
                      name: {
                        attribute_type: 'string'
                      }
                    }
                  }
                }
              }
            }
          })
        }
      );

      await waitFor(() => screen.getByTitle('Add model'));

      expect(screen.queryByText('testModel')).toBeFalsy();

      const addModelButton = screen.getByTitle('Add model');
      addModelButton.click();

      const select = screen.getByRole('button', {name: /[^Add].*/});
      fireEvent.mouseDown(select);

      const options = within(screen.getByRole('listbox'));

      const testModel = options.getByText('testModel');
      testModel.click();

      const addBtn = screen.getByRole('button', {name: 'Add'});
      addBtn.click();

      await waitFor(() => screen.queryByRole('tab'));

      expect(screen.getByText('testModel')).toBeTruthy();
    });

    it('remove model', async () => {
      render(
        <RedcapForm
          config={{
            testModel: {
              scripts: [
                {
                  attributes: {
                    name: 'redcap_field_name'
                  }
                }
              ],
              each: ['record'],
              invert: true,
              identifier_fields: ['non_standard_field']
            }
          }}
          project_name='test'
          job={{
            name: 'redcap loader',
            schema: {},
            params: {},
            secrets: {}
          }}
          update={(newConfig: Config) => {
            expect(newConfig).toEqual({});
          }}
        />,
        {
          wrapper: redcapSpecWrapper({
            props: {
              models: {
                testModel: {
                  documents: {},
                  template: {
                    attributes: {
                      name: {
                        attribute_type: 'string'
                      }
                    }
                  }
                }
              }
            }
          })
        }
      );

      await waitFor(() => screen.getByTitle('Add model'));

      expect(screen.queryByText('testModel')).toBeTruthy();

      const removeModelButton = screen.getByRole('button', {
        name: 'Remove model'
      });
      fireEvent.click(removeModelButton);
    });

    it('add script', async () => {
      const initialConfig = {
        testModel: {
          scripts: [
            {
              attributes: {
                name: 'redcap_field_name'
              }
            }
          ],
          each: ['record'],
          invert: true,
          identifier_fields: ['non_standard_field']
        }
      };

      render(
        <RedcapForm
          config={initialConfig}
          project_name='test'
          job={{
            name: 'redcap loader',
            schema: {},
            params: {},
            secrets: {}
          }}
          update={(newConfig: Config) => {
            expect(newConfig).toEqual({
              testModel: {
                ...initialConfig.testModel,
                scripts: [
                  {
                    attributes: {}
                  },
                  ...initialConfig.testModel.scripts
                ]
              }
            });
          }}
        />,
        {
          wrapper: redcapSpecWrapper({
            props: {
              models: {
                testModel: {
                  documents: {},
                  template: {
                    attributes: {
                      name: {
                        attribute_type: 'string'
                      }
                    }
                  }
                }
              }
            }
          })
        }
      );

      await waitFor(() => screen.getByTitle('Add model'));

      expect(screen.queryByText('testModel')).toBeTruthy();

      const addScriptBtn = screen.getByRole('button', {
        name: 'Add Script'
      });
      addScriptBtn.click();
    });

    it('remove script', async () => {
      const initialConfig = {
        testModel: {
          scripts: [
            {
              attributes: {
                name: 'redcap_field_name'
              }
            }
          ],
          each: ['record'],
          invert: true,
          identifier_fields: ['non_standard_field']
        }
      };

      render(
        <RedcapForm
          config={initialConfig}
          project_name='test'
          job={{
            name: 'redcap loader',
            schema: {},
            params: {},
            secrets: {}
          }}
          update={(newConfig: Config) => {
            expect(newConfig).toEqual({
              testModel: {
                ...initialConfig.testModel,
                scripts: []
              }
            });
          }}
        />,
        {
          wrapper: redcapSpecWrapper({
            props: {
              models: {
                testModel: {
                  documents: {},
                  template: {
                    attributes: {
                      name: {
                        attribute_type: 'string'
                      }
                    }
                  }
                }
              }
            }
          })
        }
      );

      await waitFor(() => screen.getByTitle('Add model'));

      expect(screen.queryByText('testModel')).toBeTruthy();

      const removeScriptBtn = screen.getByRole('button', {
        name: 'Delete script'
      });
      removeScriptBtn.click();
    });

    it('copy scripts', async () => {
      const initialConfig = {
        testModel: {
          scripts: [
            {
              attributes: {
                name: 'redcap_field_name'
              }
            }
          ],
          each: ['record'],
          invert: true,
          identifier_fields: ['non_standard_field']
        }
      };

      render(
        <RedcapForm
          config={initialConfig}
          project_name='test'
          job={{
            name: 'redcap loader',
            schema: {},
            params: {},
            secrets: {}
          }}
          update={(newConfig: Config) => {
            expect(newConfig).toEqual({
              testModel: {
                ...initialConfig.testModel,
                scripts: [
                  initialConfig.testModel.scripts[0],
                  initialConfig.testModel.scripts[0]
                ]
              }
            });
          }}
        />,
        {
          wrapper: redcapSpecWrapper({
            props: {
              models: {
                testModel: {
                  documents: {},
                  template: {
                    attributes: {
                      name: {
                        attribute_type: 'string'
                      }
                    }
                  }
                }
              }
            }
          })
        }
      );

      await waitFor(() => screen.getByTitle('Add model'));

      expect(screen.queryByText('testModel')).toBeTruthy();

      const copyScriptBtn = screen.getByRole('button', {
        name: 'Copy script'
      });
      copyScriptBtn.click();
    });

    it('set identifier_fields', async () => {
      const initialConfig = {
        testModel: {
          scripts: [
            {
              attributes: {
                name: 'redcap_field_name'
              }
            }
          ],
          each: ['record'],
          invert: true,
          identifier_fields: ['non_standard_field']
        }
      };

      render(
        <RedcapForm
          config={initialConfig}
          project_name='test'
          job={{
            name: 'redcap loader',
            schema: {},
            params: {},
            secrets: {}
          }}
          update={(newConfig: Config) => {
            expect(newConfig).toEqual({
              testModel: {
                ...initialConfig.testModel,
                identifier_fields: ['non_standard_field', 'another_field']
              }
            });
          }}
        />,
        {
          wrapper: redcapSpecWrapper({
            props: {
              models: {
                testModel: {
                  documents: {},
                  template: {
                    attributes: {
                      name: {
                        attribute_type: 'string'
                      }
                    }
                  }
                }
              }
            }
          })
        }
      );

      await waitFor(() => screen.getByTitle('Add model'));

      expect(screen.queryByText('testModel')).toBeTruthy();

      const identifierInput = screen.getByDisplayValue('non_standard_field');
      fireEvent.change(identifierInput, {
        target: {value: 'non_standard_field,another_field'}
      });
    });

    it('edit scripts', async () => {
      const initialConfig = {
        testModel: {
          scripts: [
            {
              attributes: {
                name: 'redcap_field_name'
              }
            }
          ],
          each: ['record'],
          invert: true,
          identifier_fields: ['non_standard_field']
        }
      };

      render(
        <RedcapForm
          config={initialConfig}
          project_name='test'
          job={{
            name: 'redcap loader',
            schema: {},
            params: {},
            secrets: {}
          }}
          update={(newConfig: Config) => {
            expect(newConfig).toEqual({
              testModel: {
                ...initialConfig.testModel,
                identifier_fields: ['non_standard_field', 'another_field']
              }
            });
          }}
        />,
        {
          wrapper: redcapSpecWrapper({
            props: {
              models: {
                testModel: {
                  documents: {},
                  template: {
                    attributes: {
                      name: {
                        attribute_type: 'string'
                      }
                    }
                  }
                }
              }
            }
          })
        }
      );

      await waitFor(() => screen.getByTitle('Add model'));

      expect(screen.queryByText('testModel')).toBeTruthy();

      const identifierInput = screen.getByDisplayValue('non_standard_field');
      fireEvent.change(identifierInput, {
        target: {value: 'non_standard_field,another_field'}
      });
    });

    it('edit scripts without affecting other scripts', async () => {
      // Not quite as accurate a test as I was hoping, since it doesn't
      //   catch the stale config issue I was facing with multiple
      //   scripts...but keeping it here in case it's useful.
      const initialConfig = {
        testModel: {
          scripts: [
            {
              attributes: {
                name: {
                  redcap_field: 'redcap_field_name',
                  value: 'age',
                  exists: true
                }
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: 'redcap_field_name',
                  value: 'age',
                  combine: ''
                }
              }
            }
          ],
          each: ['record'],
          invert: true,
          identifier_fields: ['non_standard_field']
        }
      };

      const state1 = {
        testModel: {
          ...initialConfig.testModel,
          scripts: [
            {
              attributes: {
                name: {
                  redcap_field: 'redcap_field_name',
                  value: 'age'
                }
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: 'redcap_field_name',
                  value: 'age',
                  combine: ''
                }
              }
            }
          ]
        }
      };

      const {rerender} = render(
        <RedcapForm
          config={initialConfig}
          project_name='test'
          job={{
            name: 'redcap loader',
            schema: {},
            params: {},
            secrets: {}
          }}
          update={(newConfig: Config) => {
            expect(newConfig).toEqual(state1);
          }}
        />,
        {
          wrapper: redcapSpecWrapper({
            props: {
              models: {
                testModel: {
                  documents: {},
                  template: {
                    attributes: {
                      name: {
                        attribute_type: 'string'
                      }
                    }
                  }
                }
              }
            }
          })
        }
      );

      await waitFor(() => screen.getByTitle('Add model'));

      expect(screen.queryByText('testModel')).toBeTruthy();

      // 3+3 remove properties
      let removePropertyBtn = screen.getAllByTitle('Remove property')[2];
      removePropertyBtn.click();

      const state2 = {
        testModel: {
          ...initialConfig.testModel,
          scripts: [
            {
              attributes: {
                name: {
                  redcap_field: 'redcap_field_name',
                  value: 'age'
                }
              }
            },
            {
              attributes: {
                name: {
                  redcap_field: 'redcap_field_name',
                  value: 'age'
                }
              }
            }
          ]
        }
      };

      rerender(
        <RedcapForm
          config={state1}
          project_name='test'
          job={{
            name: 'redcap loader',
            schema: {},
            params: {},
            secrets: {}
          }}
          update={(newConfig: Config) => {
            expect(newConfig).toEqual(state2);
          }}
        />
      );

      // 2 + 3 remove properties
      removePropertyBtn = screen.getAllByTitle('Remove property')[4];
      removePropertyBtn.click();
    });
  });
});
