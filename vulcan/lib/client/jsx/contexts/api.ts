import {MutableRefObject, useEffect, useRef, useState} from "react";
import {VulcanState} from "../reducers/vulcan_reducer";
import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
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
    postInputs(workflowName: string, session: VulcanSession): Promise<SessionStatusResponse> {
        return Promise.resolve(defaultSessionStatusResponse);
    },
    getWorkflows(): Promise<WorkflowsResponse> {
        return Promise.resolve(defaultWorkflowsResponse);
    },
}

export function useApi(): typeof defaultApiHelpers {
    const invoke = useActionInvoker();

    const [workCount, setWorkCount] = useState(0);
    const workRef = useRef(workCount);

    const vulcanPath = (endpoint: string) => `${CONFIG.vulcan_host}${endpoint}`;
    const vulcanPost = (endpoint: string, params: Object) => {
        return fetch(endpoint, {
            method: 'POST',
            credentials: 'include',
            headers: headers('json'),
            body: JSON.stringify({
                ...params
            })
        }).then(checkStatus);
    };

    const rawVulcanGet = (endpoint: string) => {
        return fetch(endpoint, {
            method: 'GET',
            credentials: 'include',
            headers: headers('json')
        });
    };

    const vulcanGet = (endpoint: string) => {
        return rawVulcanGet(endpoint).then(checkStatus);
    };

    // a
    const getWorkflows = (): Promise<WorkflowsResponse> => {
        return vulcanGet(vulcanPath(ROUTES.fetch_workflows()))
            .then(handleFetchSuccess)
            .catch(handleFetchError);
    };

    // b
    const postInputs = (workflowName: string, session: VulcanSession): Promise<SessionStatusResponse> => {
        return vulcanPost(vulcanPath(ROUTES.submit(workflowName)), session);
    }

    // c
    const getData = (url: string) => {
        return vulcanGet(url).then(handleFetchSuccess).catch(handleFetchError).then(data => {
            // TODO: In the future, we should set content type headers to inform the client, for now we aggressively
            // try to parse JSON
            try {
                return JSON.parse(data);
            } catch {
                return data;
            }
        })
    };

    return {
        isLoading: workCount > 0,
        scheduleWork<T>(work: Promise<T>): Promise<T> {
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
        },
        getData,
        getWorkflows,
        postInputs,
    }
}