import {useEffect, useCallback, MutableRefObject} from 'react';
import {SessionStatusResponse, VulcanFigure, Workflow} from '../api_types';
import {VulcanState} from '../reducers/vulcan_reducer';
import {
  allWorkflowPrimaryInputSources,
  inputValueNonEmpty
} from '../selectors/workflow_selectors';

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
    : `${project_name}/figure/new/${workflow_name}`;

export const defaultSessionStorageHelpers = {
  getLocalSession(
    workflow_name: string,
    project_name: string,
    figure_id: number
  ): Promise<SessionStatusResponse['session'] & VulcanFigure | null> {
    return Promise.resolve(null);
  }
};

// invoked every time
export function useLocalSessionStorage(
  state: VulcanState,
  props: {storage?: typeof localStorage} = {}
): typeof defaultSessionStorageHelpers {
  const storage = props.storage || localStorage;
  const {session, workflow, figure} = state;

  useEffect(() => {
    if (storage && session && figure) {
      storage.setItem(localStorageKey({
        figure_id: figure.figure_id,
        workflow_name: session.workflow_name,
        project_name: session.project_name
      }), JSON.stringify({...session, ...figure}));
    }
  }, [session, storage, figure]);

  const getLocalSession = useCallback(
    (workflow_name: string, project_name: string, figure_id: number) => {
      let storedSession: any = storage.getItem(
        localStorageKey({workflow_name, figure_id, project_name})
      );
      if (!storedSession) return Promise.resolve(null);

      try {
        const parsedSession: VulcanState['session'] & VulcanFigure = JSON.parse(storedSession);
        if (parsedSession.project_name !== project_name)
          return Promise.resolve(null);
        return Promise.resolve(parsedSession);
      } catch (e) {
        // No guarantees that the stored session is really valid, gracefully clear that session in that case.
        console.error(e);
        return Promise.resolve(null);
      }
    },
    [storage]
  );

  return {
    getLocalSession
  };
}
