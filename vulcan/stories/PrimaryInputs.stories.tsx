import React, {useContext, useEffect} from 'react';
import { ComponentStory, ComponentMeta } from '@storybook/react';

import PrimaryInputs from "../lib/client/jsx/components/workflow/session/primary_inputs";
import {VulcanContext} from "../lib/client/jsx/contexts/vulcan_context";
import {setWorkflow} from "../lib/client/jsx/actions/vulcan_actions";
import {workflowsResponse} from "../lib/client/jsx/test_utils/fixtures/workflows-response";

interface Parameterization {
}

function ParameterizedPrimaryInputs(args: Parameterization) {
  const {dispatch} = useContext(VulcanContext);

  useEffect(() => {
    dispatch(setWorkflow(
      workflowsResponse.workflows[0],
      (workflowsResponse.workflows[0].projects as any)[0]
    ))
  });

  return <PrimaryInputs/>;
}

export default {
  title: 'Inputs/PrimaryInputs',
  component: ParameterizedPrimaryInputs,
} as ComponentMeta<typeof ParameterizedPrimaryInputs>;

const Template: ComponentStory<typeof ParameterizedPrimaryInputs> = (args) => <ParameterizedPrimaryInputs {...args}/>;

export const Test = Template.bind({});
Test.args = {};
