import {useEffect, useCallback} from 'react';
import {VulcanStorage} from '../api_types';
import {VulcanState} from '../reducers/vulcan_reducer';

const localStorageKey = ({
  workspaceId,
  projectName,
}: {
  workspaceId: number;
  projectName: string;
}) =>
  `${projectName}/workspace/${workspaceId}`

export const defaultSessionStorageHelpers = {
  getLocalSession(
    workspaceId: number, projectName: string
  ): Promise<VulcanStorage | null> {
    return Promise.resolve(null);
  },
  clearLocalSession(): void {}
};

// invoked every time
export function useLocalSessionStorage(
  state: VulcanState,
  props: {storage?: typeof localStorage} = {}
): typeof defaultSessionStorageHelpers {
  const storage = props.storage || localStorage;
  const {workspace, workflow, workspaceId, status} = state;

  useEffect(() => {
    if (storage && workflow && workspace && workspaceId && status) {
      storage.setItem(
        localStorageKey({
          workspaceId: workspaceId,
          projectName: state.projectName
        }),
        // What should this look like???
        JSON.stringify({workspace: {...workspace}, status: {...status}})
      );
    }
  }, [storage, workflow, workspace, workspaceId, status]);

  const getLocalSession = useCallback(
    (workspaceId: number, projectName: string) => {
      let storedSession: any = storage.getItem(
        localStorageKey({workspaceId, projectName})
      );
      if (!storedSession) return Promise.resolve(null);

      try {
        const parsedSession: VulcanStorage = JSON.parse(storedSession);
        if (parsedSession.workspace.project !== projectName)
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

  const clearLocalSession = () => {
    if (workspaceId) {
      storage.setItem(
        localStorageKey({
          workspaceId: workspaceId,
          projectName: state.projectName
        }),
        JSON.stringify({})
      );
    };
  };

  return {
    getLocalSession,
    clearLocalSession
  };
}
