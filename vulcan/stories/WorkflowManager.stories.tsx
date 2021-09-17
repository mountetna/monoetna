import {ComponentMeta, ComponentStory} from "@storybook/react";
import React, {useContext, useEffect} from "react";
import {Workflow, WorkflowsResponse} from "../lib/client/jsx/api_types";
import {workflowUtilsBuilder} from "../lib/client/jsx/test_utils/workflow_utils";
import {MockContext} from "../lib/client/jsx/storybook-mocks";
import WorkflowManager from "../lib/client/jsx/components/workflow/workflow_manager";
import {delay} from "../lib/client/jsx/contexts/session_sync";
import {createFakeBackend, Script} from "../lib/client/jsx/test_utils/mock-backend";
import {DataEnvelope} from "../lib/client/jsx/components/workflow/user_interactions/inputs/input_types";

function MockedWorkflowManager({
  workflow,
  scripts = {},
  projectName = (workflow.projects ? workflow.projects[0] : undefined) || "",
  workflowName = workflow.name,
}: {
  workflow: Workflow,
  scripts?: DataEnvelope<Script>
  projectName?: string,
  workflowName?: string
}) {
  const {setOverrides} = useContext(MockContext);

  useEffect(() => {
    setOverrides(createFakeBackend([workflow], scripts))
  }, [scripts, setOverrides, workflow])

  return <WorkflowManager workflowName={workflowName} projectName={projectName}/>;
}

export default {
  title: 'Integration/WorkflowManager', component: MockedWorkflowManager,
} as ComponentMeta<typeof MockedWorkflowManager>;

const Template: ComponentStory<typeof MockedWorkflowManager> = (args) => <MockedWorkflowManager {...args}/>;

export const ExampleWorkflow = Template.bind({});
ExampleWorkflow.args = {
  workflow: workflowUtilsBuilder()
    .do(utils => utils.setWorkflow('test'))
    .do(utils => utils.addPrimaryInput('1__some_input', {label: "With labels, too!"}))
    .do(utils => utils.addPrimaryInput('1__other_input'))
    .do(utils => utils.addPrimaryInput('2__other_other_input'))
    .do(utils => utils.addPrimaryInput('ungrouped_input'))
    .do(utils => utils.addStep('join-strings-1', {
      in: [
        { source: '1__some_input', id: 'a' },
        { source: 'ungrouped_input', id: 'b' },
      ],
      out: ['joined'],
      run: 'scripts/join-strings',
    }))
    .do(utils => utils.addStep('join-strings-2', {
      in: [
        { source: 'join-strings-1/joined', id: 'a' },
        { source: '1__other_input', id: 'b' },
      ],
      out: ['joined'],
      run: 'scripts/join-strings',
    }))
    .do(utils => utils.addStep('checkboxes', {
      in: [
        { source: 'join-strings-1/joined', id: 'a' },
        { source: 'join-strings-2/joined', id: 'b' }
      ],
      out: ['checked'],
      run: 'ui-queries/checkboxes',
    }))
    .do(utils => utils.addStep('show_it', {
      in: [ { source: 'checkboxes/checked', id: 'a' }, { source: '2__other_other_input', id: 'b' } ],
      out: [],
      run: 'ui-outputs/raw',
    }))
    .workflow,
  scripts: {
    'scripts/join-strings': async ({a, b}: DataEnvelope<string>) => {
      await delay(4000);
      return {
        joined: a + b,
      }
    },
  },
};
