import React, {useEffect} from 'react';
import { ComponentStory, ComponentMeta } from '@storybook/react';
import PrimaryInputs from "../lib/client/jsx/components/workflow/session/primary_inputs";
import {useWorkflowUtils} from "../lib/client/jsx/test_utils/workflow_utils";
import {defaultWorkflow, WorkflowInput} from "../lib/client/jsx/api_types";

interface Parameterization {
  inputs: {[k: string]: Partial<WorkflowInput>}
}

function ParameterizedPrimaryInputs({inputs = {}}: Parameterization) {
  const utils = useWorkflowUtils();

  useEffect(() => {
    utils.setWorkflow('test');
    Object.entries(inputs).forEach(([name, options]) => {
      utils.addPrimaryInput(name, options);
    })
  }, [inputs, utils]);

  return <PrimaryInputs/>;
}

export default {
  title: 'Inputs/PrimaryInputs',
  component: ParameterizedPrimaryInputs,
} as ComponentMeta<typeof ParameterizedPrimaryInputs>;

const Template: ComponentStory<typeof ParameterizedPrimaryInputs> = (args) => <ParameterizedPrimaryInputs {...args}/>;

export const AllWithoutDefaults = Template.bind({});
AllWithoutDefaults.args = {
  inputs: {
    'int': { type: 'int' },
    'float': { type: 'float' },
    'boolean': { type: 'boolean' },
    'string': { type: 'string' },
  }
};

export const AllWithDefaults = Template.bind({});
AllWithDefaults.args = {
  inputs: {
    'int': { type: 'int', default: -22 },
    'float': { type: 'float', default: 50.3 },
    'boolean': { type: 'boolean', default: true },
    'string': { type: 'string', default: "Quack!" },
  }
};

export const WithLabelAndDocs = Template.bind({});
WithLabelAndDocs.args = {
  inputs: {
    'string': { type: 'string', label: "The String of Fate", doc: "Some sort of longer, doc string would go here." },
  }
};
