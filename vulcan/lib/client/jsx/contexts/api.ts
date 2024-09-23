import {useCallback, useRef, useState} from 'react';
import {showMessages} from 'etna-js/actions/message_actions';
import {
  checkStatus,
  handleFetchError,
  handleFetchSuccess,
  headers
} from 'etna-js/utils/fetch';
import {
  defaultSessionStatusResponse,
  defaultWorkflowsResponse,
  FiguresResponse,
  SessionStatusResponse,
  VulcanFigureSession,
  VulcanSession,
  WorkflowsResponse
} from '../api_types';

export const defaultApiHelpers = {
  showErrors<T>(work: Promise<T>): Promise<T> {
    work.catch((e) => {
      console.error(e);
    });

    return work;
  },
  getData(url: string): Promise<any> {
    return new Promise(() => null);
  },
  postInputs(session: VulcanSession): Promise<SessionStatusResponse> {
    return new Promise(() => null);
  },
  pollStatus(session: VulcanSession): Promise<SessionStatusResponse> {
    return new Promise(() => null);
  },
  getWorkflows(): Promise<WorkflowsResponse> {
    return new Promise(() => null);
  },
  fetchFigures(projectName: string): Promise<FiguresResponse> {
    return new Promise(() => null);
  },
  fetchFigure(
    projectName: string,
    figureId: number
  ): Promise<VulcanFigureSession> {
    return new Promise(() => null);
  },
  createFigure(projectName: string, params: any): Promise<VulcanFigureSession> {
    return new Promise(() => null);
  },
  updateFigure(projectName: string, params: any): Promise<VulcanFigureSession> {
    return new Promise(() => null);
  },
  deleteFigure(
    projectName: string,
    figureId: number
  ): Promise<VulcanFigureSession> {
    return new Promise(() => null);
  },
  updateFigureDependencies(projectName: string, figureId: number): Promise<VulcanFigureSession> {
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

  const vulcanDelete = useCallback((endpoint: string) => {
    return fetch(endpoint, {
      method: 'DELETE',
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

  const getWorkflows = useCallback((): Promise<WorkflowsResponse> => {
    // old: ROUTES.fetch_workflows
    return vulcanGet(vulcanPath('/api/workflows'))
      .then(handleFetchSuccess)
      .catch(handleFetchError);
  }, [vulcanGet, vulcanPath]);

  const postInputs = useCallback(
    (session: VulcanSession): Promise<SessionStatusResponse> => {
      if (!session.workflow_name) {
        return Promise.reject(
          new Error('No workflow selected, bug in client.')
        );
      }

      return vulcanPost(
        // old: ROUTES.submit
        vulcanPath(`/api/${session.project_name}/session/${workflow_name}`),
        session
      )
        .then(handleFetchSuccess)
        .catch(handleFetchError);
    },
    [vulcanPath, vulcanPost]
  );

  const pollStatus = useCallback(
    (session: VulcanSession): Promise<SessionStatusResponse> => {
      if (!session.workflow_name) {
        return Promise.reject(
          new Error('No workflow selected, bug in client.')
        );
      }

      return vulcanPost(
        // Old: ROUTES.status
        vulcanPath(`/api/${session.project_name}/session/${session.workflow_name}/status`),
        session
      )
        .then(handleFetchSuccess)
        .catch(handleFetchError);
    },
    [vulcanPath, vulcanPost]
  );

  const getData = useCallback(
    (url: string) => {
      return vulcanGet(url)
        .then(handleFetchSuccess)
        .catch(handleFetchError)
        .then((data) => {
          // TODO: In the future, we should set content type headers to inform the client, for now we aggressively
          // try to parse JSON
          try {
            return JSON.parse(data);
          } catch {
            return data;
          }
        });
    },
    [vulcanGet]
  );

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

  const fetchFigures = useCallback(
    (projectName: string): Promise<FiguresResponse> => {
      // Old: ROUTES.fetch_figures
      return vulcanGet(vulcanPath(`/api/${projectName}/figures`))
        .then(handleFetchSuccess)
        .catch(handleFetchError);
    },
    [vulcanGet, vulcanPath]
  );

  const fetchFigure = useCallback(
    (projectName: string, figureId: number): Promise<VulcanFigureSession> => {
      // Old: ROUTES.fetch_figure
      return vulcanGet(vulcanPath(`/api/${projectName}/figure/${figureId}`))
        .then(handleFetchSuccess)
        .catch(handleFetchError);
    },
    [vulcanGet, vulcanPath]
  );

  const updateFigure = useCallback(
    (projectName: string, params: any): Promise<VulcanFigureSession> => {
      return vulcanPost(
        // Old: ROUTES.update_figure
        vulcanPath(`/api/${projectName}/figure/${params.figure_id}/update`),
        params
      )
        .then(handleFetchSuccess)
        .catch(handleFetchError);
    },
    [vulcanPost, vulcanPath]
  );

  const updateFigureDependencies = useCallback(
    (projectName: string, figure_id: number): Promise<VulcanFigureSession> => {
      return vulcanPost(
        // Old: ROUTES.update_figure
        vulcanPath(`/api/${projectName}/figure/${figure_id}/update`),
        {
          update_dependencies: true,
          comment: 'Updating dependencies and snapshot'
        }
      )
        .then(handleFetchSuccess)
        .catch(handleFetchError);
    },
    [vulcanPost, vulcanPath]
  );

  const createFigure = useCallback(
    (projectName: string, params: any): Promise<VulcanFigureSession> => {
      // Old: ROUTES.create_figure
      return vulcanPost(vulcanPath(`/api/${projectName}/figure/create`), params)
        .then(handleFetchSuccess)
        .catch(handleFetchError);
    },
    [vulcanPost, vulcanPath]
  );

  const deleteFigure = useCallback(
    (projectName: string, figureId: number): Promise<VulcanFigureSession> => {
      return vulcanDelete(
        // Old: ROUTES.delete_figure
        vulcanPath(`/api/${projectName}/figure/${figureId}`)
      )
        .then(handleFetchSuccess)
        .catch(handleFetchError);
    },
    [vulcanPath, vulcanDelete]
  );

  return {
    showErrors,
    getData,
    getWorkflows,
    postInputs,
    pollStatus,
    deleteFigure,
    createFigure,
    updateFigure,
    fetchFigures,
    fetchFigure,
    updateFigureDependencies
  };
}
