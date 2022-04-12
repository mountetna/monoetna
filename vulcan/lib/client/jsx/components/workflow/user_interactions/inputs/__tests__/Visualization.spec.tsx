import React from 'react';
import NestedSelectAutocompleteInput from '../nested_select_autocomplete';
import {DataEnvelope} from '../input_types';
import {
  integrateElement,
  setupBefore
} from '../../../../../test_utils/integration';
import {Maybe, some} from '../../../../../selectors/maybe';
import {act, ReactTestInstance} from 'react-test-renderer';
import {matchesTypePredicate} from '../../../../../test_utils/rendered';
import Autocomplete from '@material-ui/lab/Autocomplete';
import { AnyPlotly, ScatterPlotly } from '../visualizations';
import StringInput from '../string';
import SelectAutocompleteInput from '../select_autocomplete';

/*
Goals:
  - [x] Offers plot type choice,
  - [x] unless given from start
  - [x] defaults filled in to 'value' when plot type chosen
  - [x] preset removes components of 'key', but sets value[key]==val
  - [ ] user input updates key'd value
  - [ ] can restrict options to continuous/discrete_cols
*/

describe('VisualizationInput', () => {
  const onChange = setupBefore(() => jest.fn());
  const value = setupBefore(() => null as Maybe<DataEnvelope<any>>);
  const data = setupBefore(() => {
    return {
      'data_frame': {
        'cont1': [1,2,3],
        'cont2': [0.1, 0.2, 0.3],
        'disc1': ["a","b","c"],
        'disc2': ["z","y","x"],
      },
      'continuous_cols': ["cont1", "cont2"],
      'discrete_cols': ["disc1", "disc2"],
      'preset': {
        'color_by': 'cont1'
      }
    } as DataEnvelope<any>;
  });
  
  const scatter_base = {
    'plot_type': 'scatter_plot',
    'x_by': null,
    'y_by': null,
    'color_by': 'make',
    'size': 5,
    'plot_title': 'make',
    'legend_title': 'make',
    'xlab': 'make',
    'ylab': 'make',
    'color_order': 'increasing',
    'order_when_continuous_color': false,
    'x_scale': 'linear',
    'y_scale': 'linear',
    'rows_use': {}
  }
  let value_when_scatter = scatter_base
  value_when_scatter.color_by = 'cont1'

  const integrated = setupBefore(() =>
    integrateElement(
      <AnyPlotly
        onChange={onChange.value}
        value={value.value}
        data={data.value}
      />
    )
  );
  
  function inputComponentExists(component: ReactTestInstance, type: any, label: string) {
    return component.findAll(matchesTypePredicate(type)).filter(
      (this_comp) => this_comp.props.label == label
    ).length > 0
  }
  
  function getInputComponent(component: ReactTestInstance, type: any, label: string) {
    return component.findAll(matchesTypePredicate(type)).filter(
      (this_comp) => this_comp.props.label == label
    )[0]
  }

  function getDropdownOptions(component: ReactTestInstance, label: string) {
    const fullInput = getInputComponent(component, SelectAutocompleteInput, label)
    return fullInput.find(matchesTypePredicate(Autocomplete)).props.options as string[];
  }

  async function clickDropdownOption(
    component: ReactTestInstance,
    label: string,
    optionIndex: number
  ) {
    const fullInput = getInputComponent(component, SelectAutocompleteInput, label)
    const autocomplete = fullInput.find(matchesTypePredicate(Autocomplete))
    
    const options = autocomplete.props.options;

    await act(async () => {
      autocomplete.props.onChange(null, options[optionIndex]);
    });
  };
  
  async function choosePlotType(component: ReactTestInstance, choice: string) {
    await clickDropdownOption(
      component,
      "Plot Type",
      getDropdownOptions(component, "Plot Type").findIndex((val) => val == choice)
    )
  }
  
  it('offers plot type choices, including the test-target-options scatter_plot', () => {
    const {node} = integrated.value;
    
    expect(onChange.value).toHaveBeenCalledWith([{"plot_type": null}])
    expect(inputComponentExists(node.root, SelectAutocompleteInput, "Plot Type")).toEqual(true)
    expect(getDropdownOptions(node.root, "Plot Type").includes('scatter_plot')).toEqual(true)
  })
  
  it('on plot choice, fills in default/preset values', async () => {
    const {node} = integrated.value;
    
    await choosePlotType(node.root, 'scatter_plot');
    
    expect(onChange.value).toHaveBeenLastCalledWith(some(value_when_scatter)) // includes many standard, and 1 preset, values.
    console.log(node.root.findAll(matchesTypePredicate(SelectAutocompleteInput)).length)
  })
  
  describe('with a plotType chosen', () => {
    value.replace(() => some(value_when_scatter));
    
    fit('renders non-preset input components', () => {
      const {node} = integrated.value;
      // Standard
      expect(inputComponentExists(node.root, SelectAutocompleteInput, 'X-Axis Data')).toEqual(true)
      // Preset
      expect(inputComponentExists(node.root, SelectAutocompleteInput, 'Color Data')).toEqual(false)
    })
  });
  
  describe('with a setPlotType', () => {
    integrated.replace(() => integrateElement(
      <ScatterPlotly
        onChange={onChange.value}
        value={value.value}
        data={data.value}
      />
    ));

    it('skips plot_type selection', () => {
      const {node} = integrated.value;

      expect(inputComponentExists(node.root, SelectAutocompleteInput, "plot_type")).toEqual(false)
      let stored_value = value.value as DataEnvelope<any>
      expect(stored_value['plot_type']).toEqual('scatter_plot')
    });
  });
  
});
