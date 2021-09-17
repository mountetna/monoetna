import {useEffect, useCallback, MutableRefObject} from "react";
import {SessionStatusResponse, Workflow} from "../api_types";
import {VulcanState} from "../reducers/vulcan_reducer";
import {allWorkflowPrimaryInputSources, inputValueNonEmpty} from "../selectors/workflow_selectors";

const localStorageKey = (workflow: Workflow, projectName: string) => `${projectName}/${workflow.name}.session`;

export const defaultSessionStorageHelpers = {
  getLocalSession(workflow: Workflow, projectName: string): Promise<SessionStatusResponse['session'] | null> {
    return Promise.resolve(null);
  }
}

export function useLocalSessionStorage(state: VulcanState,
  props: { storage?: typeof localStorage } = {}
): typeof defaultSessionStorageHelpers {

  const storage = props.storage || localStorage;
  const {session, workflow} = state;

  useEffect(() => {
    if (workflow && workflow.name && storage && session.workflow_name === workflow.name) {
      storage.setItem(localStorageKey(workflow, session.project_name), JSON.stringify(session));
    }
  }, [session, workflow, storage])

  const getLocalSession = useCallback((workflow: Workflow, projectName: string) => {
    let storedSession: any = storage.getItem(localStorageKey(workflow, projectName));
    if (!storedSession) return Promise.resolve(null);

    try {
      const parsedSession: VulcanState['session'] = JSON.parse(storedSession);
      if (parsedSession.project_name !== projectName) return Promise.resolve(null);
      return Promise.resolve(parsedSession);
    } catch (e) {
      // No guarantees that the stored session is really valid, gracefully clear that session in that case.
      console.error(e);
      return Promise.resolve(null);
    }
  }, [storage]);

  return {
    getLocalSession,
  }
}