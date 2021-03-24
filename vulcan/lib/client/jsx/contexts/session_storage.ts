import {useEffect, useCallback, MutableRefObject} from "react";
import {SessionStatusResponse, Workflow} from "../api_types";
import {VulcanState} from "../reducers/vulcan_reducer";

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
        if (workflow && workflow.name)
            storage.setItem(
                localStorageKey(workflow),
                JSON.stringify(session)
            );
    }, [state.session, state.workflow])

    const getLocalSession = useCallback((workflow: Workflow) => {
        let storedSession: any = storage.getItem(localStorageKey(workflow));
        if (!storedSession) return Promise.resolve(null);

        // We now need to check if the input names have changed.
        // If all the workflow's primary inputs are NOT present
        //   in the stored session, we'll return `null` and get
        //   a new session.
        storedSession = JSON.parse(storedSession);
        if (!allInputsDefined(workflow, storedSession.inputs)) {
            storage.removeItem(localStorageKey(workflow));
            return Promise.resolve(null);
        }

        return Promise.resolve(storedSession);
    }, [storage]);

    return {
        getLocalSession,
    }
}