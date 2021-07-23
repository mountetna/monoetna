import {setDownloadedData, setStatus, setWorkflow, VulcanAction} from "../actions/vulcan_actions";
import {
  defaultStepStatus,
  defaultWorkflow,
  defaultWorkflowInput,
  defaultWorkflowStep,
  StatusString, StepStatus,
  Workflow,
  WorkflowInput,
  WorkflowStep
} from "../api_types";
import {createStepStatusFixture, createUpdatedStatusFixture} from "./fixtures";
import VulcanReducer, {defaultVulcanState, VulcanState} from "../reducers/vulcan_reducer";
import {VulcanContext} from "../contexts/vulcan_context";
import {useContext, useState} from "react";
import {splitSource, statusOfStep} from "../selectors/workflow_selectors";
import {Step} from "@material-ui/core";

export function useWorkflowUtils(): WorkflowUtils {
  const {dispatch, stateRef} = useContext(VulcanContext);
  const [utils] = useState(() => {
    return new WorkflowUtils(dispatch, stateRef);
  });
  return utils;
}

export function workflowUtilsBuilder() {
  const stateRef = { current: defaultVulcanState };
  return new WorkflowUtils(dispatch, stateRef);

  function dispatch(action: VulcanAction) {
    stateRef.current = VulcanReducer(stateRef.current, action);
  }
}

export class WorkflowUtils {
  public workflow: Readonly<Workflow> = defaultWorkflow;
  public steps: { [k: string]: WorkflowStep } = {};

  constructor(private dispatch: (action: VulcanAction) => void, private stateRef: { current: VulcanState },) {
  }

  get status() {
    return this.stateRef.current.status;
  }

  // Useful for chaining
  // workflowUtils.do(utils => ).do(utils => );
  do(f: (utils: WorkflowUtils) => void): WorkflowUtils {
    f(this);
    return this;
  }

  setWorkflow(name: string,
    workflow: Partial<Workflow> = {},
    projects = workflow.projects || this.workflow.projects || [name],
    projectName = projects[0] || "test"
  ) {
    this.workflow = {...this.workflow, ...workflow, name, projects};
    this.dispatch(setWorkflow(this.workflow, projectName));
  }

  addStep(name: string, attributes: Partial<WorkflowStep> = {}) {
    const step = {...defaultWorkflowStep, name, ...attributes};
    ({name} = step);

    let workflow = this.workflow;

    if (!(name in this.steps)) {
      workflow = {...workflow, steps: [[...workflow.steps[0], step]]};
    }

    this.steps[name] = step;
    workflow = {...workflow, steps: [workflow.steps[0].map(s => s.name !== name ? s : step)]};
    this.setWorkflow(workflow.name, workflow);
    return step;
  }

  setStatus(step: WorkflowStep | string, statusPart: StatusString | Partial<StepStatus>) {
    const stepName = typeof step !== "string" ? step.name : step;
    const stepStatus: Partial<StepStatus> = typeof statusPart === "string" ? {status: statusPart} : statusPart;
    const status = createStepStatusFixture({name: stepName, ...stepStatus});

    const updatedStatus = createUpdatedStatusFixture(this.workflow,
      this.stateRef.current.status,
      status
    );
    this.dispatch(setStatus(updatedStatus, false));

    return status;
  }

  addPrimaryInput(name: string, options: Partial<WorkflowInput> = {}) {
    let workflow = {...this.workflow, inputs: {...this.workflow.inputs}};
    let input = {...defaultWorkflowInput, ...options};
    workflow.inputs[name] = input;
    this.setWorkflow(workflow.name, workflow);
    return input;
  }


  /*
     Note -- this method will force a setStatus update to a 'complete' status with either the existing or a new
     url constructed for the purpose of this sourceName.
   */
  forceDownloadedData(sourceName: string | [string, string], value: any) {
    const workflow = this.workflow;
    let url: string | undefined;
    const [stepName, outputName] = typeof sourceName === "string" ? splitSource(sourceName) : sourceName;
    if (!stepName) throw new Error('Cannot setData for primary inputs');

    let status = statusOfStep(stepName, this.stateRef.current.status) || defaultStepStatus;
    const {downloads} = status;
    if (downloads) {
      url = downloads[outputName]
    }
    if (!url) url = "https://" + sourceName;

    // TODO: Set the inputs hash
    this.dispatch(setStatus(createUpdatedStatusFixture(workflow, this.stateRef.current.status, createStepStatusFixture({
      ...status, name: stepName, status: 'complete', downloads: {...status.downloads, [outputName]: url}
    })), false));

    this.dispatch(setDownloadedData(url, value));
  }
}