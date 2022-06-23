import {VulcanAction} from '../actions/vulcan_actions';
import {
  defaultSessionStatusResponse,
  SessionStatusResponse,
  StepStatus,
  Workflow,
  WorkflowsResponse,
  VulcanFigure
} from '../api_types';
import {
  allExpectedOutputSources,
  filterEmptyValues,
  selectFigure,
  selectSession,
  statusOfStep,
  stepOfStatus
} from '../selectors/workflow_selectors';
import {mapSome, Maybe, some, withDefault} from '../selectors/maybe';

export type DownloadedData = any; // TODO: improve typing here.
export type DownloadedStepDataMap = {[k: string]: DownloadedData};

const defaultWorkflows: WorkflowsResponse['workflows'] = [];
const defaultWorkflow: Workflow | null = null;
const defaultStatus: SessionStatusResponse['status'] = [[]];
const defaultData: DownloadedStepDataMap = {};
const defaultInputs: SessionStatusResponse['session']['inputs'] = {};
export const defaultSession: SessionStatusResponse['session'] = {
  project_name: '',
  workflow_name: '',
  key: '',
  inputs: {},
  reference_figure_id: null
};
const defaultFigure: VulcanFigure = {
  figure_id: null,
  inputs: {},
  id: null
};
const defaultValidationErrors: [string | null, string, string[]][] = [];

export const defaultVulcanState = {
  workflows: defaultWorkflows,
  workflow: defaultWorkflow as Workflow | null,
  status: defaultStatus,
  data: defaultData,

  // A subset of all steps that have buffered changes.
  bufferedSteps: [] as (string | null)[],
  // If a buffered step was filled in while remaining 'pending'
  committedStepPending: false,

  // Step marked for auto-passing by user
  autoPassSteps: [] as (string | null)[],
  triggerRun: [] as (string | null)[],

  session: defaultSession,
  outputs: defaultSessionStatusResponse.outputs,
  validationErrors: defaultValidationErrors,
  pollingState: 0,
  figure: defaultFigure
};

export type VulcanState = Readonly<typeof defaultVulcanState>;

export default function VulcanReducer(
  state: VulcanState,
  action: VulcanAction
): VulcanState {
  state = state || defaultVulcanState;

  switch (action.type) {
    case 'SET_WORKFLOWS':
      return {
        ...state,
        workflows: action.workflows
      };
    case 'MODIFY_POLLING':
      return {
        ...state,
        pollingState: Math.max(state.pollingState + action.delta, 0)
      };
    case 'SET_WORKFLOW':
      const workflowProjects = action.workflow.projects;
      if ( workflowProjects===undefined || (workflowProjects !== null && !workflowProjects.includes(action.projectName)) )
        return state;

      return {
        ...state,
        workflow: action.workflow,
        data: defaultData
      };
    case 'SET_STATUS':
      // When a submitting step is given, filter stale inputs that result from submitting a change
      // to that step's outputs in session.inputs.
      state = filterStaleInputs(state, action);
      return {...state, status: action.status};

    case 'SET_BUFFERED_INPUT':
      if (state.bufferedSteps.includes(action.step)) {
        return state;
      }
      return {...state, bufferedSteps: [...state.bufferedSteps, action.step]};

    case 'CLEAR_BUFFERED_INPUT':
      const idx = state.bufferedSteps.indexOf(action.step);
      if (idx === -1) return state;
      const bufferedSteps = [...state.bufferedSteps];
      bufferedSteps.splice(idx, 1);
      return {...state, bufferedSteps};

    case 'SET_AUTO_PASS_STEP':
      if (state.autoPassSteps.includes(action.step)) {
        return state;
      }
      return {...state, autoPassSteps: [...state.autoPassSteps, action.step]};

    case 'CLEAR_AUTO_PASS_STEP':
      const index = state.autoPassSteps.indexOf(action.step);
      if (index === -1) return state;
      const autoPassSteps = [...state.autoPassSteps];
      autoPassSteps.splice(index, 1);
      return {...state, autoPassSteps};

    case 'SET_RUN_TRIGGER':
      if (state.triggerRun.includes(action.step)) {
        return state;
      }
      return {...state, triggerRun: [...state.triggerRun, action.step]};

    case 'CLEAR_RUN_TRIGGERS':
      const triggerRun = [...state.triggerRun].filter((val) => !action.steps.includes(val));
      return {...state, triggerRun};

    case 'SET_DOWNLOAD':
      return {
        ...state,
        data: {
          ...state.data,
          [action.url]: action.data
        }
      };

    case 'SET_SESSION':
      // Ignore sessions not for this workflow, project combination. safety guard.
      // The project name and workflow names are locked in via the set workflow action.
      const sessionWorkflow = state.workflow;
      if (
        !sessionWorkflow ||
        action.session.workflow_name !== sessionWorkflow.name
      ) {
        console.warn(
          'Cannot set session, project name / workflow name does not match.',
          sessionWorkflow,
          action
        );
        return state;
      }

      return {
        ...state,
        session: {...state.session, ...action.session},
        bufferedSteps: []
      };

    case 'REMOVE_DOWNLOADS':
      return {
        ...state,
        status: [
          state.status[0].map((status) =>
            action.stepNames.includes(status.name)
              ? {...status, downloads: null}
              : status
          )
        ]
      };

    case 'REMOVE_INPUTS':
      if (state.pollingState) {
        console.error('cannot change inputs while polling...');
        return state;
      }

      const removedInputs = {...state.session.inputs};
      action.inputs.forEach((source) => delete removedInputs[source]);

      return {
        ...state,
        session: {...state.session, inputs: filterEmptyValues(removedInputs)}
      };

    case 'SET_INPUTS':
      if (state.pollingState) {
        console.error('cannot change inputs while polling...');
        return state;
      }

      return {
        ...state,
        session: {...state.session, inputs: filterEmptyValues(action.inputs)}
      };

    case 'ADD_VALIDATION_ERRORS':
      return {
        ...state,
        validationErrors: [
          ...state.validationErrors,
          [action.stepName, action.inputLabel, action.errors]
        ]
      };

    case 'REMOVE_VALIDATION_ERRORS':
      return {
        ...state,
        validationErrors: state.validationErrors.filter(
          ([_1, _2, e]) => e !== action.errors
        )
      };

    case 'CHECK_CHANGES_READY':
      const stepName = withDefault(action.step, null);
      let ready: boolean;
      if (stepName == null) {
        ready = false;
      } else {
        const stepNum = state.status[0].findIndex(
          (element) => element['name'] == stepName
        );
        ready = state.status[0][stepNum]['status'] == 'pending';
      }
      return {
        ...state,
        committedStepPending: ready
      };

    case 'CLEAR_CHANGES_READY':
      return {
        ...state,
        committedStepPending: false
      };

    case 'SET_SESSION_AND_FIGURE':
      const currentWorkflow = state.workflow;
      if (
        !currentWorkflow ||
        action.session.workflow_name !== currentWorkflow.name
      ) {
        console.warn(
          'Cannot set session, project name / workflow name does not match.',
          currentWorkflow,
          action
        );
        return state;
      }
      return {
        ...state,
        session: {...action.session},
        figure: {...action.figure},
        bufferedSteps: []
      };

    default:
      return state;
  }
}

function filterStaleInputs(
  state: VulcanState,
  action: {type: 'SET_STATUS'} & {
    status: [StepStatus[]];
    submittingStep: Maybe<string>;
  }
) {
  const {workflow} = state;
  if (!workflow) return state;

  return withDefault(
    mapSome(action.submittingStep, (submittingStep) => {
      let newState = {...state};
      let newInputs = state.session.inputs;
      let newData = state.data;

      const hashesOfSteps = {} as {[k: string]: string};
      action.status[0].forEach((stepStatus) => {
        hashesOfSteps[stepStatus.name] = stepStatus.hash;
      });

      for (let stepStatus of state.status[0]) {
        // The submitting step is pushing a new value from the client up, thus
        // it should not have its input made stale.
        if (stepStatus.name === submittingStep) continue;

        if (hashesOfSteps[stepStatus.name] !== stepStatus.hash) {
          const step = stepOfStatus(stepStatus, workflow);
          if (!step) continue;
          allExpectedOutputSources(step).forEach((outputSource) => {
            if (newState === state) {
              newState = {
                ...state,
                session: {...state.session, inputs: {...state.session.inputs}}
              };
              newInputs = state.session.inputs;
            }
            delete newInputs[outputSource];
          });

          if (stepStatus.downloads) {
            Object.values(stepStatus.downloads).forEach((url) => {
              if (newState === state) {
                newState = {...state, data: {...state.data}};
                newData = newState.data;
              }
              delete newData[url];
            });
          }
        }
      }

      return newState;
    }),
    state
  );
}
