import {VulcanState} from "../reducers/vulcan_reducer";
import {useCallback} from "react";
import {defaultSessionSyncHelpers} from "./session_sync";

export const defaultInputStateManagement = {
  commitSessionInputChanges() {
    return Promise.resolve(true);
  },

  startInputChange(source: string) {
    return Promise.resolve(true);
  },

  cancelInputChange() {
  },

  onInputChange(source: string, value: any) {
  }
}

export function useInputStateManagement(
  stateRef: { current: VulcanState },
  requestPoll: typeof defaultSessionSyncHelpers.requestPoll,
): typeof defaultInputStateManagement {
  const commitSessionInputChanges = useCallback(() => {
    if (stateRef.current.curEditingInput == null) {
      return Promise.resolve(false);
    }

    const preCommitStepStatuses = stateRef.current.status[0]
    requestPoll();
  }, [requestPoll, stateRef]);

  return {}
}
