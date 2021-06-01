import {useCallback, useRef, useState} from "react";
import {showMessages} from "etna-js/actions/message_actions";
import {checkStatus, handleFetchError, handleFetchSuccess, headers} from "etna-js/utils/fetch";
import {
    defaultSessionStatusResponse,
    defaultWorkflowsResponse,
    SessionStatusResponse,
    VulcanSession,
    WorkflowsResponse
} from "../api_types";

export const defaultApiHelpers = {
    isLoading: false,
    scheduleWork<T>(work: Promise<T>): Promise<T> {
        work.catch((e) => {
            console.error(e);
        })

        return work;
    },
    getData(url: string): Promise<any> {
        return Promise.resolve({});
    },
    postInputs(session: VulcanSession): Promise<SessionStatusResponse> {
        return new Promise((resolve, reject) => null);
    },
    pollStatus(session: VulcanSession): Promise<SessionStatusResponse> {
        return new Promise((resolve, reject) => null);
    },
    getWorkflows(): Promise<WorkflowsResponse> {
        return new Promise((resolve, reject) => null);
    },
}

export function useApi(invoke: (a: {type: string}) => any): typeof defaultApiHelpers {
    const [workCount, setWorkCount] = useState(0);
    const workRef = useRef(workCount);

    const vulcanPath = useCallback((endpoint: string) => `${CONFIG.vulcan_host}${endpoint}`, []);
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

    const vulcanGet = useCallback((endpoint: string) => {
        return rawVulcanGet(endpoint).then(checkStatus);
    }, [rawVulcanGet]);

    const getWorkflows = useCallback((): Promise<WorkflowsResponse> => {
        return vulcanGet(vulcanPath(ROUTES.fetch_workflows()))
            .then(handleFetchSuccess)
            .catch(handleFetchError);
    }, [vulcanGet, vulcanPath]);

    const postInputs = useCallback((session: VulcanSession): Promise<SessionStatusResponse> => {
        if (!session.workflow_name) {
            return Promise.reject(new Error("No workflow selected, bug in client."));
        }

        return vulcanPost(vulcanPath(ROUTES.submit(session.project_name, session.workflow_name)), session);
    }, [vulcanPath, vulcanPost]);

    const pollStatus = useCallback((session: VulcanSession): Promise<SessionStatusResponse> => {
        if (!session.workflow_name) {
            return Promise.reject(new Error("No workflow selected, bug in client."));
        }

        return vulcanPost(vulcanPath(ROUTES.status(session.project_name, session.workflow_name)), session);
    }, [vulcanPath, vulcanPost])

    const getData = useCallback((url: string) => {
        return vulcanGet(url).then(handleFetchSuccess).catch(handleFetchError).then(data => {
            // TODO: In the future, we should set content type headers to inform the client, for now we aggressively
            // try to parse JSON
            try {
                return JSON.parse(data);
            } catch {
                return data;
            }
        })
    }, [vulcanGet]);

    const scheduleWork = useCallback(<T>(work: Promise<T>): Promise<T> => {
        setWorkCount(++workRef.current);

        work.catch(e => {
            if (!(e instanceof Array)) {
                e = [`${e}`];
            }

            console.error(e);
            invoke(showMessages(e));
        }).finally(() => {
            setWorkCount(--workRef.current);
        });

        return work;
    }, [invoke]);

    return {
        isLoading: workCount > 0,
        scheduleWork,
        getData,
        getWorkflows,
        postInputs,
        pollStatus,
    }
}