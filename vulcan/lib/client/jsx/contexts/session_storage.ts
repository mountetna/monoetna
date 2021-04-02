import {useEffect, useCallback, MutableRefObject} from "react";
import {SessionStatusResponse, Workflow} from "../api_types";
import {VulcanState} from "../reducers/vulcan_reducer";
import {allWorkflowPrimaryInputSources, inputValueNonEmpty} from "../selectors/workflow_selectors";

const localStorageKey = (workflow: Workflow) => `${workflow.name}.session`;

export const defaultSessionStorageHelpers = {
    getLocalSession(workflow: Workflow): Promise<SessionStatusResponse['session'] | null> {
        return Promise.resolve(null);
    }
}

export function useLocalSessionStorage(
    state: VulcanState,
    props: {storage?: typeof localStorage} = {}): typeof defaultSessionStorageHelpers {

    const storage = props.storage || localStorage;

    useEffect(() => {
        const {workflow, session} = state;
        if (workflow && workflow.name && storage && session.workflow_name === workflow.name) {
            storage.setItem(
                localStorageKey(workflow),
                JSON.stringify(session)
            );
        }
    }, [state.session, state.workflow])

    const getLocalSession = useCallback((workflow: Workflow) => {
        let storedSession: any = storage.getItem(localStorageKey(workflow));
        if (!storedSession) return Promise.resolve(null);

        try {
            const parsedSession: VulcanState['session'] = JSON.parse(storedSession);
            return Promise.resolve(parsedSession);
        } catch(e) {
            // No guarantees that the stored session is really valid, gracefully clear that session in that case.
            console.error(e);
            return Promise.resolve(null);
        }
    }, []);

    return {
        getLocalSession,
    }
}