import {showMessages, dismissMessages} from 'etna-js/actions/message_actions';
import {
  checkStatus,
  handleFetchError,
  handleFetchSuccess,
  headers,
  parseJSON
} from 'etna-js/utils/fetch';
import {
  AccountingReturn,
  FileContentResponse,
  MultiFileContentResponse,
  MultiFileContent,
  FlatParams,
  RunReturn,
  RunStatus,
  WorkflowsResponse,
  Workspace,
  WorkspaceStatus,
  WorkspacesResponse,
  WorkspaceRaw,
  WorkflowCreateResponse,
  CreateWorkspaceResponse,
  isRunningReturn,
  WorkspacesResponseRaw,
  LatencyReturn,
  StateReturn,
  ConfigReturn
} from '../api_types';
import { paramValuesToRaw, workspacesFromResponse } from '../selectors/workflow_selectors';
import { isSome } from '../selectors/maybe';
import { val_wrap } from '../components/workspace/ui_definitions/inputs/pieces/user_input_pieces';
import * as _ from 'lodash';

const vulcanPath = (endpoint: string) => `${CONFIG.vulcan_host}${endpoint}`;

const vulcanPostRaw = (endpoint: string, params: Object) => {
  return fetch(endpoint, {
    method: 'POST',
    credentials: 'include',
    headers: headers('json'),
    body: JSON.stringify({
      ...params
    })
  })
};

const vulcanPost = (endpoint: string, params: Object) => {
  return vulcanPostRaw(endpoint, params)
  .then(checkStatus)
  .then(handleFetchSuccess)
  .catch(handleFetchError);
};

const rawVulcanGet = (endpoint: string) => {
  return fetch(endpoint, {
    method: 'GET',
    credentials: 'include',
    headers: headers('json')
  });
};

const vulcanGet = (endpoint: string) => {
  return rawVulcanGet(endpoint)
    .then(checkStatus)
    .then(handleFetchSuccess)
    .catch(handleFetchError);
};

const vulcanDelete = (endpoint: string) => {
  return fetch(endpoint, {
    method: 'DELETE',
    credentials: 'include',
    headers: headers('json')
  });
};

const showError = (e: any, dismissOld: boolean = false) => {
  if (dismissOld) invoke(dismissMessages());
  if (!(e instanceof Array)) {
    e = [`${e}`];
  }
  console.error(e);
  showMessages(e);
};

const showErrors = <T>(work: Promise<T>, additional: (e: any) => void = (e) => {}): Promise<T> => {
  work.catch((e) => {
    if (!(e instanceof Array)) {
      e = [`${e}`];
    }

    console.error(e);
    invoke(showMessages(e));
    additional(e);
  });

  return work;
};

const getWorkflows = (projectName: string): Promise<WorkflowsResponse> => {
  return vulcanGet(vulcanPath(`/api/v2/${projectName}/workflows`))
};

const getDag = (projectName: string, workspaceId: string): Promise<WorkflowsResponse> => {
  return vulcanGet(vulcanPath(`/api/v2/${projectName}/workspace/${workspaceId}/dag`))
};

const createWorkflow = (projectName: string, repoUrl: string, workflowName: string): Promise<WorkflowCreateResponse> => {
  return vulcanPost(
    vulcanPath(`/api/v2/${projectName}/workflows/create`),
    {
      project_name: projectName,
      repo_url: repoUrl,
      workflow_name: workflowName
    });
};

const createWorkspace = (projectName: string, workflowId: number, workspaceName: string, git_request: string): Promise<CreateWorkspaceResponse> => {
  return vulcanPost(
    vulcanPath(`/api/v2/${projectName}/workspace/create`),
    {
      workflow_id: workflowId,
      git_request: git_request,
      workspace_name: workspaceName,
    });
};

const getWorkspaces = (projectName: string): Promise<WorkspacesResponse> => {
    return vulcanGet(vulcanPath(`/api/v2/${projectName}/workspace`))
      .then((val: unknown) => workspacesFromResponse(val as WorkspacesResponseRaw));
};

const getWorkspace = (projectName: string, workspaceId: number): Promise<WorkspaceRaw> => {
    // Old: ROUTES.fetch_figure
    return vulcanGet(vulcanPath(`/api/v2/${projectName}/workspace/${workspaceId}`));
};

const updateWorkspace = (projectName: string, workspaceId: number, name?: string, tags?: string[]): Promise<WorkspaceRaw> => {
  const params = {};
  if (!!name) params['name'] = name;
  if (!!tags) params['tags'] = tags;
  if (Object.keys(params).length < 1) {
    showError('UI Error: updateWorkspace was called without any updates to send.')
  }
  return vulcanPost(
    vulcanPath(`/api/v2/${projectName}/workspace/${workspaceId}/update`),
    params
  );
};

const deleteWorkspace = (projectName: string, workspaceId: number): Promise<Response> => {
    return vulcanDelete(vulcanPath(`/api/v2/${projectName}/workspace/${workspaceId}`));
};

const getFileNames = (projectName: string, workspaceId: number): Promise<{files: string[]}> => {
    return vulcanGet(vulcanPath(`/api/v2/${projectName}/workspace/${workspaceId}/file`))
};

const writeFiles = (projectName: string, workspaceId: number, files_content: MultiFileContent) => {
  return vulcanPost(
    vulcanPath(`/api/v2/${projectName}/workspace/${workspaceId}/file/write`),
    {files:
      Object.entries(files_content).map(([name, contents]) => {
        return {filename: name, content: JSON.stringify(contents)}
      })
    }
  );
};

const readFiles = (projectName: string, workspaceId: number, fileNames: string[]): Promise<MultiFileContentResponse> => {
  return vulcanPost(
    vulcanPath(`/api/v2/${projectName}/workspace/${workspaceId}/file/read`),
    {
      file_names: fileNames
    }
  );
};

const getImage = (projectName: string, workspaceId: number, imageFile: string): Promise<Response> => {
  return vulcanPostRaw(
    vulcanPath(`/api/v2/${projectName}/workspace/${workspaceId}/image/read`),
    {
      file_name: imageFile
    }
  );
};

const setConfig = (projectName: string, workspaceId: number, params: FlatParams, uiFilesSent: string[], paramsChanged: string[]): Promise<AccountingReturn> => {
    return vulcanPost(
      vulcanPath(`/api/v2/${projectName}/workspace/${workspaceId}/config`),
      {
        params,
        uiFilesSent,
        paramsChanged
      }
    );
};

const postUIValues = (projectName: string, workspaceId: number, status: WorkspaceStatus, step: string): Promise<AccountingReturn> => {
  // Only ever per a single 'step'.
  // Can assume the values have already been validated.
  // If a "inputUI", send targeted file outputs
  let uiFilesSent: string[] = [];
  let paramsChanged: string[] = [];
  let paramsUse: FlatParams = status.last_params;
  if (step in status.ui_contents) {
    let filesContent: MultiFileContent = {};
    Object.entries(status.ui_contents[step]).map(([key, val]) => {
      if (isSome(val)) {
        // These should never really not be a 'some' / at least [null] as will have been deliberately set to a (possibly null) value.
        uiFilesSent.push(key);
        filesContent[key] = val[0];
      } else {
        return Promise.reject(
          new Error(`${key} had no value, so cannot send file contents.`)
        );
      }
    })
    return showErrors(writeFiles(
      projectName,
      workspaceId,
      {...filesContent}
    ))
    .then(() => {
      return showErrors(setConfig(
        projectName,
        workspaceId,
        paramsUse,
        uiFilesSent,
        paramsChanged
      ))
    })
  } else {
    paramsUse = paramValuesToRaw(status.params)
    paramsChanged = Object.keys(status.params[step]).filter((name) => !_.isEqual(paramsUse[name], status.last_params[name]))
    return showErrors(setConfig(
      projectName,
      workspaceId,
      paramsUse,
      uiFilesSent,
      paramsChanged
    ));
  }
};

const getState = (projectName: string, configId: number): Promise<StateReturn> => {
    return vulcanGet(
      vulcanPath(`/api/v2/${projectName}/config/${configId}/state`)
    );
};

const getConfig = (projectName: string, workspaceId: number, configId: number): Promise<ConfigReturn> => {
    return vulcanGet(
      vulcanPath(`/api/v2/${projectName}/workspace/${workspaceId}/config/${configId}`)
    );
};

const requestRun = (projectName: string, workspaceId: number, configId: number): Promise<RunReturn> => {
    return vulcanPost(
      vulcanPath(`/api/v2/${projectName}/workspace/${workspaceId}/run/${configId}`), {}
    );
};

const getIsRunning = (projectName: string, workspaceId: number): Promise<isRunningReturn> => {
    return vulcanGet(vulcanPath(`/api/v2/${projectName}/workspace/${workspaceId}/running`))
};

const pullRunStatus = (projectName: string, workspaceId: number, runId: number): Promise<RunStatus> => {
    return vulcanGet(vulcanPath(`/api/v2/${projectName}/workspace/${workspaceId}/run/${runId}`))
};

const getConnectionLatency = (): Promise<LatencyReturn> => {
  return vulcanGet(vulcanPath(`/api/v2/cluster-latency`))
};

export {
  vulcanPath,
  showError,
  showErrors,
  createWorkflow,
  getWorkflows,
  createWorkspace,
  getWorkspaces,
  getWorkspace,
  updateWorkspace,
  deleteWorkspace,
  getFileNames,
  getDag,
  writeFiles,
  readFiles,
  setConfig,
  postUIValues,
  getState,
  getConfig,
  requestRun,
  getIsRunning,
  pullRunStatus,
  getImage,
  getConnectionLatency
};
