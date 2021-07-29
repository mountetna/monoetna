import React, {useContext, useEffect} from 'react';
import {ComponentMeta, ComponentStory} from "@storybook/react";
import StepOutput from "../lib/client/jsx/components/workflow/steps/step_output";
import {Workflow} from "../lib/client/jsx/api_types";
import {stepOfStatus, uiOutputOfStep} from "../lib/client/jsx/selectors/workflow_selectors";
import {VulcanContext} from "../lib/client/jsx/contexts/vulcan_context";
import {setWorkflow} from "../lib/client/jsx/actions/vulcan_actions";
import {useWorkflowUtils, workflowUtilsBuilder} from "../lib/client/jsx/test_utils/workflow_utils";

function MockStepOutput({
  workflow,
  stepName = (workflow.steps[0].find(step => !!uiOutputOfStep(step))?.name || ""),
  projectName = (workflow.projects ? workflow.projects[0] : undefined) || "",
}: {
  workflow: Workflow,
  stepName?: string,
  projectName?: string,
}) {
  const step = stepOfStatus(stepName, workflow);
  if (!step) throw new Error(`Step ${stepName} does not exist on the given workflow`);

  const utils = useWorkflowUtils();
  const {dispatch} = useContext(VulcanContext);

  useEffect(() => {
    utils.setWorkflow(workflow.name, workflow, [projectName]);
    utils.forceDownloadedData('load-data/a', require('./umap.plotly.json'));
    utils.setStatus('plot', 'complete');
  }, [dispatch, projectName, utils, workflow]);

  return <StepOutput step={step}/>
}

export default {
  title: 'Outputs/StepOutput', component: MockStepOutput,
} as ComponentMeta<typeof MockStepOutput>;

const Template: ComponentStory<typeof MockStepOutput> = (args) => <MockStepOutput {...args}/>;

export const PlotWithEffects = Template.bind({});
PlotWithEffects.args = {
  workflow: workflowUtilsBuilder()
    .do(utils => utils.addStep('load-data', {run: 'scripts/load', out: ['a']}))
    .do(utils => utils.addStep('plot', {run: 'ui-outputs/plotly', in: [{source: 'load-data/a', id: 'a'}]}))
    .workflow
};