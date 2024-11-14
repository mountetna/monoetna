import {useCallback, useRef, useState} from 'react';
import {showMessages} from 'etna-js/actions/message_actions';
import {
  checkStatus,
  handleFetchError,
  handleFetchSuccess,
  headers
} from 'etna-js/utils/fetch';
import {
  AccountingReturn,
  ConfigContent,
  defaultSessionStatusResponse,
  defaultWorkflowsResponse,
  FiguresResponse,
  FileContent,
  MultiFileContent,
  ParamsContent,
  RunReturn,
  RunStatus,
  SessionStatusResponse,
  VulcanFigureSession,
  VulcanSession,
  WorkflowsResponse,
  WorkspaceResponse,
  WorkspacesResponse
} from '../api_types';

export const defaultApiHelpers = {
  showErrors<T>(work: Promise<T>): Promise<T> {
    work.catch((e) => {
      console.error(e);
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
  createWorkflow(projectName: string, repoUrl: string, workflowName: string): Promise<Response> {
    return new Promise(() => null);
  },
  getWorkflows(projectName: string): Promise<Response | WorkflowsResponse> {
    return new Promise(() => null);
  },
  createWorkspace(projectName: string, workflowId: number, branch: string, workspaceName: string): Promise<Response | WorkspaceResponse> {
    return new Promise(() => null);
  },
  getWorkspaces(projectName: string): Promise<Response | WorkspacesResponse> {
    return new Promise(() => null);
  },
  getWorkspace(projectName: string, workspaceId: number): Promise<Response | WorkspaceResponse> {
    return new Promise(() => null);
  },
  getFileNames(projectName: string, workspaceId: number): Promise<Response | string[]> {
    return new Promise(() => null);
  },
  writeFiles(projectName: string, workspaceId: number, content: FileContent | MultiFileContent) {
    return new Promise(() => null);
  },
  readFiles(projectName: string, workspaceId: number, fileNames: string[]): Promise<Response | MultiFileContent> {
    return new Promise(() => null);
  },
  setConfig(projectName: string, workspaceId: number, params: ParamsContent): Promise<Response | AccountingReturn> {
    return new Promise(() => null);
  },
  requestRun(projectName: string, workspaceId: number, configId: number): Promise<Response | RunReturn> {
    return new Promise(() => null);
  },
  pullRunStatus(projectName: string, workspaceId: number, runId: number): Promise<Response | RunStatus> {
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
  const vulcanPost = useCallback((endpoint: string, params: Object) => {
    return fetch(endpoint, {
      method: 'POST',
      credentials: 'include',
      headers: headers('json'),
      body: JSON.stringify({
        ...params
      })
    }).then(checkStatus);
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

  const rawVulcanGetWithParams = useCallback((endpoint: string, params: Object) => {
    return fetch(endpoint, {
      method: 'GET',
      credentials: 'include',
      headers: headers('json'),
      body: JSON.stringify({
        ...params
      })
    });
  }, []);

  const vulcanGetWithParams = useCallback(
    (endpoint: string, params: Object) => {
      return rawVulcanGetWithParams(endpoint, params)
        .then(checkStatus)
        .then(handleFetchSuccess)
        .catch(handleFetchError);
    },
    [rawVulcanGetWithParams]
  );

  const vulcanDelete = useCallback((endpoint: string) => {
    return fetch(endpoint, {
      method: 'DELETE',
      credentials: 'include',
      headers: headers('json')
    });
  }, []);

  const getWorkflows = useCallback((projectName: string): Promise<Response | WorkflowsResponse> => {
    // old: ROUTES.fetch_workflows
    return vulcanGet(vulcanPath(`/api/v2/${projectName}/workflows`))
  }, [vulcanGet, vulcanPath]);

  const createWorkflow = useCallback((projectName: string, repoUrl: string, workflowName: string): Promise<Response> => {
    return vulcanPost(
      vulcanPath(`/api/v2/${projectName}/workflows/create`),
      {
        project_name: projectName,
        repo_url: repoUrl,
        workflow_name: workflowName
      }).then(handleFetchSuccess).catch(handleFetchError);
  }, [vulcanPost, vulcanPath]);

  // const postInputs = useCallback(
  //   (session: VulcanSession): Promise<SessionStatusResponse> => {
  //     if (!session.workflow_name) {
  //       return Promise.reject(
  //         new Error('No workflow selected, bug in client.')
  //       );
  //     }

  //     return vulcanPost(
  //       // old: ROUTES.submit
  //       vulcanPath(`/api/${session.project_name}/session/${workflow_name}`),
  //       session
  //     )
  //       .then(handleFetchSuccess)
  //       .catch(handleFetchError);
  //   },
  //   [vulcanPath, vulcanPost]
  // );

  // const pollStatus = useCallback(
  //   (session: VulcanSession): Promise<SessionStatusResponse> => {
  //     if (!session.workflow_name) {
  //       return Promise.reject(
  //         new Error('No workflow selected, bug in client.')
  //       );
  //     }

  //     return vulcanPost(
  //       // Old: ROUTES.status
  //       vulcanPath(`/api/${session.project_name}/session/${session.workflow_name}/status`),
  //       session
  //     )
  //       .then(handleFetchSuccess)
  //       .catch(handleFetchError);
  //   },
  //   [vulcanPath, vulcanPost]
  // );

  // const getData = useCallback(
  //   (url: string) => {
  //     return vulcanGet(url)
  //       .then(handleFetchSuccess)
  //       .catch(handleFetchError)
  //       .then((data) => {
  //         // TODO: In the future, we should set content type headers to inform the client, for now we aggressively
  //         // try to parse JSON
  //         try {
  //           return JSON.parse(data);
  //         } catch {
  //           return data;
  //         }
  //       });
  //   },
  //   [vulcanGet]
  // );

  const createWorkspace = useCallback((projectName: string, workflowId: number, branch: string, workspaceName: string): Promise<Response | WorkspaceResponse> => {
    return vulcanPost(
      vulcanPath(`/api/v2/${projectName}/workspace/create`),
      {
        workflow_id: workflowId,
        branch: branch,
        workspace_name: workspaceName,
        git_version: "v1",
      }).then(handleFetchSuccess).catch(handleFetchError);
  }, [vulcanPost, vulcanPath]);

  const showErrors = useCallback(
    <T>(work: Promise<T>): Promise<T> => {
      work.catch((e) => {
        if (!(e instanceof Array)) {
          e = [`${e}`];
        }

        console.error(e);
        invoke(showMessages(e));
      });

      return work;
    },
    [invoke]
  );

  const getWorkspaces = useCallback(
    (projectName: string): Promise<Response | WorkspacesResponse> => {
      // Old: ROUTES.fetch_figures
      return vulcanGet(vulcanPath(`/api/v2/${projectName}/workspace`));
  }, [vulcanGet, vulcanPath]);

  const getWorkspace = useCallback(
    (projectName: string, workspaceId: number): Promise<Response | WorkspaceResponse> => {
      // Old: ROUTES.fetch_figure
      return vulcanGet(vulcanPath(`/api/v2/${projectName}/workspace/${workspaceId}`));
  }, [vulcanGet, vulcanPath]);

  const getFileNames = useCallback(
    (projectName: string, workspaceId: number): Promise<Response | string[]> => {
      return vulcanGet(vulcanPath(`/api/v2/${projectName}/workspace/${workspaceId}/file`))
  }, [vulcanGet, vulcanPath]);

  const writeFiles = useCallback(
    (projectName: string, workspaceId: number, content: FileContent | MultiFileContent) => {
      content = Array.isArray(content) ? content : [content]
      return vulcanPost(
        vulcanPath(`/api/v2/${projectName}/workspace/${workspaceId}/file/write`),
        {
          files: content
        });
  }, [vulcanPost, vulcanPath]);

  const readFiles = useCallback(
    (projectName: string, workspaceId: number, fileNames: string[]): Promise<Response| MultiFileContent> => {
      return vulcanGetWithParams(
        vulcanPath(`/api/v2/${projectName}/workspace/${workspaceId}/file/read`),
        {
          file_names: fileNames
        });
  }, [vulcanGetWithParams, vulcanPath]);

  const setConfig = useCallback(
    (projectName: string, workspaceId: number, params: ParamsContent): Promise<Response | AccountingReturn> => {
      return vulcanPost(
        vulcanPath(`/api/v2/${projectName}/workspace/${workspaceId}/config`),
        params
      );
  }, [vulcanPost, vulcanPath]);

  const requestRun = useCallback(
    (projectName: string, workspaceId: number, configId: number): Promise<Response | RunReturn> => {
      return vulcanPost(
        vulcanPath(`/api/v2/${projectName}/workspace/${workspaceId}/run/${configId}`), {}
      );
  }, [vulcanPost, vulcanPath]);

  const pullRunStatus = useCallback(
    (projectName: string, workspaceId: number, runId: number): Promise<Response | RunStatus> => {
      return vulcanGet(vulcanPath(`/api/v2/${projectName}/workspace/${workspaceId}/run/${runId}`))
  }, [vulcanGet, vulcanPath]);

  // const updateFigure = useCallback(
  //   (projectName: string, params: any): Promise<VulcanFigureSession> => {
  //     return vulcanPost(
  //       // Old: ROUTES.update_figure
  //       vulcanPath(`/api/${projectName}/figure/${params.figure_id}/update`),
  //       params
  //     )
  //       .then(handleFetchSuccess)
  //       .catch(handleFetchError);
  //   },
  //   [vulcanPost, vulcanPath]
  // );

  // const updateFigureDependencies = useCallback(
  //   (projectName: string, figure_id: number): Promise<VulcanFigureSession> => {
  //     return vulcanPost(
  //       // Old: ROUTES.update_figure
  //       vulcanPath(`/api/${projectName}/figure/${figure_id}/update`),
  //       {
  //         update_dependencies: true,
  //         comment: 'Updating dependencies and snapshot'
  //       }
  //     )
  //       .then(handleFetchSuccess)
  //       .catch(handleFetchError);
  //   },
  //   [vulcanPost, vulcanPath]
  // );

  // Not yet implemented
  // const deleteFigure = useCallback(
  //   (projectName: string, figureId: number): Promise<VulcanFigureSession> => {
  //     return vulcanDelete(
  //       // Old: ROUTES.delete_figure
  //       vulcanPath(`/api/${projectName}/figure/${figureId}`)
  //     )
  //       .then(handleFetchSuccess)
  //       .catch(handleFetchError);
  //   },
  //   [vulcanPath, vulcanDelete]
  // );

  return {
    showErrors,
    //getData,
    createWorkflow,
    getWorkflows,
    // postInputs,
    // pollStatus,
    createWorkspace,
    getWorkspaces,
    getWorkspace,
    getFileNames,
    writeFiles,
    readFiles,
    setConfig,
    requestRun,
    pullRunStatus
  };
}
