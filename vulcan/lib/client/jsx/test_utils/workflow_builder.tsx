import {Dispatch} from "react";
import {setStatus, setWorkflow, VulcanAction} from "../actions/vulcan";
import {defaultWorkflow, defaultWorkflowStep, StatusString, StepStatus, Workflow, WorkflowStep} from "../api_types";
import {createStatusFixture, createStepStatusFixture, createUpdatedStatusFixture} from "./fixtures";
import {VulcanState} from "../reducers/vulcan_reducer";
import {VulcanContextData} from "../contexts/vulcan_context";

export class WorkflowBuilder {
  public workflow: Readonly<Workflow> = defaultWorkflow;
  public steps: { [k: string]: WorkflowStep } = {};

  constructor(private dispatch: (action: VulcanAction) => Promise<void>, private stateRef: { current: VulcanState },) {
  }

  // combine with the result of integrateElement.
  static fromSetup({dispatch, contextData}: { dispatch: (action: VulcanAction) => Promise<void>, contextData: VulcanContextData }) {
    return new WorkflowBuilder(dispatch, contextData.stateRef);
  }

  async setWorkflow(name: string, workflow: Partial<Workflow> = {},
    projects = workflow.projects || this.workflow.projects || [],
    projectName = projects[0] || "test"
  ) {
    this.workflow = {...this.workflow, ...workflow, name, projects: workflow.projects || [name]};
    await this.dispatch(setWorkflow(this.workflow, projectName));
  }

  async addStep(name: string, attributes: Partial<WorkflowStep> = {}) {
    const step = {...defaultWorkflowStep, name, ...attributes};
    ({name} = step);

    let workflow = this.workflow;

    if (!(name in this.steps)) {
      workflow = {...workflow, steps: [[...workflow.steps[0], step]]};
    }

    this.steps[name] = step;
    workflow = {...workflow, steps: [workflow.steps[0].map(s => s.name !== name ? s : step)] };
    await this.setWorkflow(workflow.name, workflow);
    return step;
  }

  async setStatus(step: WorkflowStep | string, status: StatusString) {
    const stepName = typeof step !== "string" ? step.name : step;
    await this.dispatch(setStatus(createUpdatedStatusFixture(
      this.workflow,
      this.stateRef.current.status,
      createStepStatusFixture({name: stepName, status})
    )));
  }
}