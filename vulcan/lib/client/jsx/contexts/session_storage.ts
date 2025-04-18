import {useEffect, useCallback} from 'react';
import {
  SessionStatusResponse,
  VulcanFigure,
  VulcanFigureSession,
  VulcanSession
} from '../api_types';
import {VulcanState} from '../reducers/vulcan_reducer';
import {cwlName} from '../selectors/workflow_selectors';

const localStorageKey = ({
  figure_id,
  project_name,
  workflow_name
}: {
  figure_id?: number | null;
  project_name: string;
  workflow_name: string;
}) =>
  figure_id
    ? `${project_name}/figure/${figure_id}`
    : `${project_name}/figure/new/${cwlName(workflow_name)}`;

export const defaultSessionStorageHelpers = {
  getLocalSession(
    workflow_name: string,
    project_name: string,
    figure_id: number | null
  ): Promise<(SessionStatusResponse['session'] & VulcanFigure) | null> {
    return Promise.resolve(null);
  },
  clearLocalSession(
    workflow_name: string,
    project_name: string,
    figure_id: number | null
  ): void {}
};

// invoked every time
export function useLocalSessionStorage(
  state: VulcanState,
  props: {storage?: typeof localStorage} = {}
): typeof defaultSessionStorageHelpers {
  const storage = props.storage || localStorage;
  const {session, workflow, figure} = state;

  useEffect(() => {
    if (storage && session && figure && workflow) {
      storage.setItem(
        localStorageKey({
          figure_id: figure.figure_id,
          workflow_name: session.workflow_name,
          project_name: session.project_name
        }),
        JSON.stringify({...figure, ...session})
      );
    }
  }, [session, storage, figure, workflow]);

  const getLocalSession = useCallback(
    (workflow_name: string, project_name: string, figure_id: number | null) => {
      let storedSession: any = storage.getItem(
        localStorageKey({workflow_name, figure_id, project_name})
      );
      if (!storedSession) return Promise.resolve(null);

      try {
        const parsedSession: VulcanFigureSession = JSON.parse(storedSession);
        if (parsedSession.project_name !== project_name)
          {return Promise.resolve(null);}
        return Promise.resolve(parsedSession);
      } catch (e) {
        // No guarantees that the stored session is really valid, gracefully clear that session in that case.
        console.error(e);
        return Promise.resolve(null);
      }
    },
    [storage]
  );

  const clearLocalSession = (
    workflow_name: string,
    project_name: string,
    figure_id: number | null
  ) => {
    storage.setItem(
      localStorageKey({
        figure_id: figure.figure_id,
        workflow_name: session.workflow_name,
        project_name: session.project_name
      }),
      JSON.stringify({})
    );
  };

  return {
    getLocalSession,
    clearLocalSession
  };
}
