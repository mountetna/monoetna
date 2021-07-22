import {ComponentMeta, ComponentStory} from "@storybook/react";
import React, {useContext, useEffect} from "react";
import {SessionStatusResponse, VulcanSession, Workflow, WorkflowsResponse} from "../lib/client/jsx/api_types";
import {WorkflowUtils, workflowUtilsBuilder} from "../lib/client/jsx/test_utils/workflow_utils";
import {MockContext} from "../lib/client/jsx/storybook-mocks";
import WorkflowManager from "../lib/client/jsx/components/workflow/workflow_manager";
import {delay} from "../lib/client/jsx/contexts/session_sync";

function MockedWorkflowManager({
  workflow,
  projectName = (workflow.projects ? workflow.projects[0] : undefined) || "",
  workflowName = workflow.name,

}: {
  workflow: Workflow,
  projectName?: string,
  workflowName?: string
}) {
  const {setOverrides} = useContext(MockContext);

  useEffect(() => {
    setOverrides({
      async getWorkflows(): Promise<WorkflowsResponse> {
        // simulate delays for more realistic feels
        await delay(1000)
        return {
          workflows: [workflow],
        }
      },

      // async pollStatus(session: VulcanSession): Promise<SessionStatusResponse> {
      //
      // }
    })
  }, [setOverrides, workflow])

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
    .do(utils => utils.addPrimaryInput('1__some_input'))
    .do(utils => utils.addPrimaryInput('1__other_input'))
    .do(utils => utils.addPrimaryInput('2__other_other_input'))
    .do(utils => utils.addPrimaryInput('ungrouped_input'))
    .do(utils => utils.addStep('do_work', {
      in: [ { source: '1__some_input', id: 'a' }, { source: 'ungrouped_input', id: 'b' } ],
      out: ['a'],
      run: 'scripts/do_work',
    }))
    .do(utils => utils.addStep('do_work_2', {
      in: [ { source: 'do_work/a', id: 'a' }, { source: '1__other_input', id: 'b' } ],
      out: ['b'],
      run: 'scripts/do_work_2',
    }))
    .do(utils => utils.addStep('query_it', {
      in: [ { source: 'do_work/a', id: 'a' }, { source: '1__other_input', id: 'b' } ],
      out: ['a'],
      run: 'ui-queries/checkboxes',
    }))
    .do(utils => utils.addStep('show_it', {
      in: [ { source: 'query_it/a', id: 'a' }, { source: '2__other_other_input', id: 'b' } ],
      out: [],
      run: 'ui-outputs/raw',
    }))
    .workflow,
};
