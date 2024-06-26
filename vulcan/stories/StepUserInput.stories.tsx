import React, {useContext, useEffect, useMemo, useState} from 'react';
import { ComponentStory, ComponentMeta } from '@storybook/react';
import PrimaryInputs from "../lib/client/jsx/components/workflow/session/primary_inputs";
import {useWorkflowUtils, workflowUtilsBuilder} from "../lib/client/jsx/test_utils/workflow_utils";
import {defaultWorkflow, TYPE, Workflow, WorkflowInput, WorkflowStep} from "../lib/client/jsx/api_types";
import { DataEnvelope } from '../lib/client/jsx/components/workflow/user_interactions/inputs/input_types';
import StepUserInput from '../lib/client/jsx/components/workflow/steps/step_user_input';
import { useWorkflow } from '../lib/client/jsx/contexts/workflow_context';
import { VulcanContext } from '../lib/client/jsx/contexts/vulcan_context';
import {mapSome, Maybe, some, withDefault} from "../lib/client/jsx/selectors/maybe";
import {isPendingUiQuery} from "../lib/client/jsx/selectors/workflow_selectors";

interface Parameterization {
  cwlParams?: DataEnvelope<any>,
  type: string
}


function ParameterizedStepUserInput({cwlParams = {}, type}: Parameterization) {
  const utils = useWorkflowUtils();
  const [step, setState] = useState(null as Maybe<WorkflowStep>);
  const {state} = useContext(VulcanContext);
  const {status, data, session} = state;

  useEffect(() => {
    utils.setWorkflow('test');
    Object.keys(cwlParams).forEach(paramName => {
      utils.addStep(paramName, { out: ['output'] });
    });
    setState(some(utils.addStep('test-inputs', {
        run: `ui-queries/${type}`,
        label: "Title here",
        doc: "doc string would show here",
        in: Object.keys(cwlParams).map(paramName => ({source: `${paramName}/output`, id: paramName})),
        out: ['result'],
    })));

    Object.keys(cwlParams).forEach(paramName => {
      utils.forceDownloadedData(`${paramName}/output`, cwlParams[paramName]);
    });
  }, [utils, cwlParams, type]);

  return withDefault(mapSome(step, step => {
    if (!isPendingUiQuery(step, status, data, session)) {
      return <div>
        Not all data inputs are loaded.
      </div>;
    }

    return <StepUserInput step={step} hideLabel={false}/>;
  }), null);
}

export default {
  title: 'Inputs/StepUserInput',
  component: ParameterizedStepUserInput,
} as ComponentMeta<typeof ParameterizedStepUserInput>;

const Template: ComponentStory<typeof ParameterizedStepUserInput> = (args: any) => <ParameterizedStepUserInput {...args}/>;

export const DiffExpSC = Template.bind({});
DiffExpSC.args = {
  type: TYPE.DIFF_EXP_SC,
  cwlParams: {
      'data_frame': require('./mockDF.json')
  }
};

export const Visualization = Template.bind({});
Visualization.args = {
  type: TYPE.ANY_VIZ,
  cwlParams: {
      'data_frame': require('./mockDF.json'),
      'continuous_cols': require('./mockDF_cont_cols.json'),
      'discrete_cols': require('./mockDF_disc_cols.json')
  }
};

export const ScatterPlotlyFULL = Template.bind({});
ScatterPlotlyFULL.args = {
  type: TYPE.SCATTER_PLOTLY,
  cwlParams: {
      'data_frame': require('./mockDF.json'),
      'continuous_cols': require('./mockDF_cont_cols.json'),
      'discrete_cols': require('./mockDF_disc_cols.json')
  }
};

export const BarPlotlyFULL = Template.bind({});
BarPlotlyFULL.args = {
  type: TYPE.BAR_PLOTLY,
  cwlParams: {
      'data_frame': require('./mockDF.json'),
      'continuous_cols': require('./mockDF_cont_cols.json'),
      'discrete_cols': require('./mockDF_disc_cols.json')
  }
};

export const YPlotlyFULL = Template.bind({});
YPlotlyFULL.args = {
  type: TYPE.Y_PLOTLY,
  cwlParams: {
      'data_frame': require('./mockDF.json'),
      'continuous_cols': require('./mockDF_cont_cols.json'),
      'discrete_cols': require('./mockDF_disc_cols.json')
  }
};

export const ScatterPlotlyUMAP = Template.bind({});
ScatterPlotlyUMAP.args = {
  type: TYPE.SCATTER_PLOTLY,
  cwlParams: {
      'data_frame': require('./mockDF.json'),
      'preset': {
        'x_by': '0', 'y_by': '1', 'color_by': 'leiden',
        'xlab': 'UMAP_1', 'ylab': 'UMAP_2',
        'hover_data': 'hover_data'}
  }
};

export const RecordSelection = Template.bind({});
RecordSelection.args = {
  type: TYPE.MULTIPLE_MULTISELECT_STRING_ALL,
  cwlParams: {
    'selection_options': {
        'Experiment': ["exp1", "expA", "exp42"],
        'Tissue': ["blood", "tumor", "dLN"],
        'Fraction': ["myeloid", "lymphoid", "CD45neg"]
      }
  }
};

export const RecordConfirmation = Template.bind({});
RecordConfirmation.args = {
  type: TYPE.CHECKBOXES,
  cwlParams: {
    'a': ['rec____________1', 'rec____________2', 'rec____________3', 'rec____________4', 'rec____________5']
  }
};

export const ColorSelection = Template.bind({});
ColorSelection.args = {
  type: TYPE.NESTED_SELECT_AUTOCOMPLETE,
  cwlParams: {
    'color_options': require('./color_options.json')
  }
};

export const BatchSelection = Template.bind({});
BatchSelection.args = {
  type: TYPE.SELECT_AUTOCOMPLETE,
  cwlParams:  {
    'batch_options': ['1', '2', '3'],
    'recommendation': ['1', '2', '3']
  }
}

export const manyOptionDropdown = Template.bind({});
manyOptionDropdown.args = {
  type: TYPE.SELECT_AUTOCOMPLETE,
  cwlParams:  {
    'options': require('./human_genes_list.json')
  }
}

export const dropdownCheckboxes = Template.bind({});
dropdownCheckboxes.args = {
  type: TYPE.SINGLE_DROPDOWN_MULTICHECKBOX,
  cwlParams:  {
    'a': {experiment: ['1', '2'], tissue: ['a', 'b']}
  }
}
    
export const DataTransformation = Template.bind({});
DataTransformation.args = {
  type: TYPE.DATA_TRANSFORMATION,
  cwlParams: {
    'data_frame': {
      'record_name_01': {
        '0': 1,
        '1': 0.25,
        '2': 'data'
      },
      'record_name_02': {
        '0': 200,
        '1': 1.111,
        '2': '=IF(A2>100, 1, 0)'
      }
    }
  }
}

export const nestedTest = Template.bind({});
nestedTest.args = {
  type: TYPE.NESTED_SELECT_AUTOCOMPLETE,
  cwlParams: {
    'options-a': {
      option1: {
        suboption1: null,
        suboption2: {
          grandchild1: null,
          grandchild2: null
        }
      }
    },
    'options-b': {
      option2: {
        another1: {
          stepchild1: null,
          stepchild2: null
        }
      }
    }
  }
};