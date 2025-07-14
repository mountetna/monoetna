import {useCallback, useRef, useState} from 'react';
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
  ClusterStatusReturn
} from '../api_types';
import { paramValuesToRaw, workspacesFromResponse } from '../selectors/workflow_selectors';
import { isSome } from '../selectors/maybe';
import { val_wrap } from '../components/workspace/ui_definitions/inputs/pieces/user_input_pieces';
import * as _ from 'lodash';

export const defaultApiHelpers = {
  vulcanPath(endpoint: string): string {
    return `${CONFIG.vulcan_host}${endpoint}`
  },
  showError(e: string, dismissOld = false) {
    console.error(e);
  },
  showErrors<T>(work: Promise<T>, additional: (e: any) => void = (e) => {}): Promise<T> {
    work.catch((e) => {
      console.error(e);
      additional(e);
    });

    return work;
  },
  // getData(url: string): Promise<any> {
  //   return new Promise(() => null);
  // },
  // postInputs(session: VulcanSession): Promise<SessionStatusResponse> {
  //   return new Promise(() => null);
  // },
  // pollStatus(session: VulcanSession): Promise<SessionStatusResponse> {
  //   return new Promise(() => null);
  // },
  createWorkflow(projectName: string, repoUrl: string, workflowName: string): Promise<WorkflowCreateResponse> {
    return new Promise(() => null);
  },
  getWorkflows(projectName: string): Promise<WorkflowsResponse> {
    return new Promise(() => null);
  },
  createWorkspace(projectName: string, workflowId: number, workspaceName: string, git_request: string): Promise<CreateWorkspaceResponse> {
    return new Promise(() => null);
  },
  getWorkspaces(projectName: string): Promise<WorkspacesResponse> {
    return new Promise(() => null);
  },
  getWorkspace(projectName: string, workspaceId: number): Promise<WorkspaceRaw> {
    return new Promise(() => null);
  },
  updateWorkspace(projectName: string, workspaceId: number, name?: string, tags?: string[]): Promise<WorkspaceRaw> {
    return new Promise(() => null);
  },
  deleteWorkspace(projectName: string, workspaceId: number): Promise<Response> {
    return new Promise(() => null);
  },
  getFileNames(projectName: string, workspaceId: number): Promise<{files: string[]}> {
    return new Promise(() => null);
  },
  writeFiles(projectName: string, workspaceId: number, content: MultiFileContent): Promise<Response> {
    return new Promise(() => null);
  },
  readFiles(projectName: string, workspaceId: number, fileNames: string[]): Promise<MultiFileContentResponse> {
    return new Promise(() => null);
  },
  getImage(projectName: string, workspaceId: number, imageFile: string): Promise<Response> {
    return new Promise(() => null);
  },
  setConfig(projectName: string, workspaceId: number, params: FlatParams, uiFilesSent: string[], paramsChanged: string[]): Promise<AccountingReturn> {
    return new Promise(() => null);
  },
  postUIValues(projectName: string, workspaceId: number, status: WorkspaceStatus, steps: string | null): Promise<AccountingReturn> {
    return new Promise(() => null);
  },
  requestRun(projectName: string, workspaceId: number, configId: number): Promise<RunReturn> {
    return new Promise(() => null);
  },
  getIsRunning(projectName: string, workspaceId: number): Promise<isRunningReturn> {
    return new Promise(() => null);
  },
  pullRunStatus(projectName: string, workspaceId: number, runId: number): Promise<RunStatus> {
    return new Promise(() => null);
  },
  getClusterLatency(projectName: string): Promise<LatencyReturn> {
    return new Promise(() => null);
  },
  getClusterStatus(projectName: string): Promise<ClusterStatusReturn> {
    return new Promise(() => null);
  }
};

export function useApi(
  invoke: (a: {type: string}) => any
): typeof defaultApiHelpers {
  const vulcanPath = useCallback(
    (endpoint: string) => `${CONFIG.vulcan_host}${endpoint}`,
    []
  );

  const vulcanPostRaw = useCallback((endpoint: string, params: Object) => {
    return fetch(endpoint, {
      method: 'POST',
      credentials: 'include',
      headers: headers('json'),
      body: JSON.stringify({
        ...params
      })
    })
  }, []);

  const vulcanPost = useCallback((endpoint: string, params: Object) => {
    return vulcanPostRaw(endpoint, params)
    .then(checkStatus)
    .then(handleFetchSuccess)
    .catch(handleFetchError);
  }, []);

  const rawVulcanGet = useCallback((endpoint: string) => {
    return fetch(endpoint, {
      method: 'GET',
      credentials: 'include',
      headers: headers('json')
    });
  }, []);

  const vulcanGet = useCallback(
    (endpoint: string) => {
      return rawVulcanGet(endpoint)
        .then(checkStatus)
        .then(handleFetchSuccess)
        .catch(handleFetchError);
    },
    [rawVulcanGet]
  );

  const vulcanDelete = useCallback((endpoint: string) => {
    return fetch(endpoint, {
      method: 'DELETE',
      credentials: 'include',
      headers: headers('json')
    });
  }, []);

  const showError = useCallback((e: any, dismissOld: boolean = false) => {
    if (dismissOld) invoke(dismissMessages());
    console.error(e);
    invoke(showMessages([e]));
  }, [invoke]);
  const showErrors = useCallback(
    <T>(work: Promise<T>, additional: (e: any) => void = (e) => {}): Promise<T> => {
      work.catch((e) => {
        if (!(e instanceof Array)) {
          e = [`${e}`];
        }

        console.error(e);
        invoke(showMessages(e));
        additional(e);
      });

      return work;
    },
    [invoke]
  );

  const getWorkflows = useCallback((projectName: string): Promise<WorkflowsResponse> => {
    return vulcanGet(vulcanPath(`/api/v2/${projectName}/workflows`))
  }, [vulcanGet, vulcanPath]);

  const createWorkflow = useCallback((projectName: string, repoUrl: string, workflowName: string): Promise<WorkflowCreateResponse> => {
    return vulcanPost(
      vulcanPath(`/api/v2/${projectName}/workflows/create`),
      {
        project_name: projectName,
        repo_url: repoUrl,
        workflow_name: workflowName
      });
  }, [vulcanPost, vulcanPath]);

  const createWorkspace = useCallback((projectName: string, workflowId: number, workspaceName: string, git_request: string): Promise<CreateWorkspaceResponse> => {
    return vulcanPost(
      vulcanPath(`/api/v2/${projectName}/workspace/create`),
      {
        workflow_id: workflowId,
        git_request: git_request,
        workspace_name: workspaceName,
      });
  }, [vulcanPost, vulcanPath]);

  const getWorkspaces = useCallback(
    (projectName: string): Promise<WorkspacesResponse> => {
      return vulcanGet(vulcanPath(`/api/v2/${projectName}/workspace`))
        .then((val: unknown) => workspacesFromResponse(val as WorkspacesResponseRaw));
  }, [vulcanGet, vulcanPath]);

  const getWorkspace = useCallback(
    (projectName: string, workspaceId: number): Promise<WorkspaceRaw> => {
      // Old: ROUTES.fetch_figure
      return vulcanGet(vulcanPath(`/api/v2/${projectName}/workspace/${workspaceId}`));
  }, [vulcanGet, vulcanPath]);

  const updateWorkspace = useCallback(
    (projectName: string, workspaceId: number, name?: string, tags?: string[]): Promise<WorkspaceRaw> => {
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
  }, [vulcanPost, vulcanPath]);

  const deleteWorkspace = useCallback(
    (projectName: string, workspaceId: number): Promise<Response> => {
      return vulcanGet(vulcanPath(`/api/v2/${projectName}/workspace/${workspaceId}/delete`));
  }, [vulcanGet, vulcanPath]);

  const getFileNames = useCallback(
    (projectName: string, workspaceId: number): Promise<{files: string[]}> => {
      return vulcanGet(vulcanPath(`/api/v2/${projectName}/workspace/${workspaceId}/file`))
  }, [vulcanGet, vulcanPath]);

  const writeFiles = useCallback(
    (projectName: string, workspaceId: number, files_content: MultiFileContent) => {
      return vulcanPost(
        vulcanPath(`/api/v2/${projectName}/workspace/${workspaceId}/file/write`),
        {files:
          Object.entries(files_content).map(([name, contents]) => {
            return {filename: name, content: JSON.stringify(contents)}
          })
        }
      );
  }, [vulcanPost, vulcanPath]);

  const readFiles = useCallback(
    (projectName: string, workspaceId: number, fileNames: string[]): Promise<MultiFileContentResponse> => {
      return vulcanPost(
        vulcanPath(`/api/v2/${projectName}/workspace/${workspaceId}/file/read`),
        {
          file_names: fileNames
        });
  }, [vulcanPost, vulcanPath]);

  const getImage = useCallback(
    (projectName: string, workspaceId: number, imageFile: string): Promise<Response> => {
      return vulcanPostRaw(
        vulcanPath(`/api/v2/${projectName}/workspace/${workspaceId}/image/read`),
        {
          file_name: imageFile
        });
  }, [vulcanPost, vulcanPath]);

  const setConfig = useCallback(
    (projectName: string, workspaceId: number, params: FlatParams, uiFilesSent: string[], paramsChanged: string[]): Promise<AccountingReturn> => {
      return vulcanPost(
        vulcanPath(`/api/v2/${projectName}/workspace/${workspaceId}/config`),
        {
          params,
          uiFilesSent,
          paramsChanged
        }
      );
  }, [vulcanPost, vulcanPath]);

  const postUIValues = useCallback(
    (projectName: string, workspaceId: number, status: WorkspaceStatus, step: string): Promise<AccountingReturn> => {
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
    },
    [vulcanPost, vulcanPath, writeFiles, setConfig]
  );

  const requestRun = useCallback(
    (projectName: string, workspaceId: number, configId: number): Promise<RunReturn> => {
      return vulcanPost(
        vulcanPath(`/api/v2/${projectName}/workspace/${workspaceId}/run/${configId}`), {}
      );
  }, [vulcanPost, vulcanPath]);

  const getIsRunning = useCallback(
    (projectName: string, workspaceId: number): Promise<isRunningReturn> => {
      return vulcanGet(vulcanPath(`/api/v2/${projectName}/workspace/${workspaceId}/running`))
  }, [vulcanGet, vulcanPath]);

  const pullRunStatus = useCallback(
    (projectName: string, workspaceId: number, runId: number): Promise<RunStatus> => {
      return vulcanGet(vulcanPath(`/api/v2/${projectName}/workspace/${workspaceId}/run/${runId}`))
  }, [vulcanGet, vulcanPath]);

  const getClusterLatency = useCallback((projectName: string): Promise<LatencyReturn> => {
    return vulcanGet(vulcanPath(`/api/v2/${projectName}/cluster-latency`))
  }, [vulcanGet, vulcanPath]);

  const getClusterStatus = useCallback((projectName: string): Promise<ClusterStatusReturn> => {
    return vulcanGet(vulcanPath(`/api/v2/${projectName}/cluster-status`))
  }, [vulcanGet, vulcanPath]);

  return {
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
    writeFiles,
    readFiles,
    setConfig,
    postUIValues,
    requestRun,
    getIsRunning,
    pullRunStatus,
    getImage,
    getClusterLatency,
    getClusterStatus,
  };
}
