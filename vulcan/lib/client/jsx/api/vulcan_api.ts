import * as _ from 'lodash';

import {
    checkStatus,
    handleFetchSuccess,
    handleFetchError,
    headers,
    isJSON
} from 'etna-js/utils/fetch';
import {shouldDownloadStepData} from '../utils/workflow';
import {SessionStatusResponse, VulcanSession, WorkflowsResponse} from "./types";
import {VulcanContext} from "../contexts/vulcan";
import {VulcanState} from "../reducers/vulcan";

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

export const getWorkflows = (): Promise<WorkflowsResponse> => {
    return vulcanGet(vulcanPath(ROUTES.fetch_workflows()))
        .then(handleFetchSuccess)
        .catch(handleFetchError);
};

export const postInputs = (workflowName: string, session: VulcanSession): Promise<SessionStatusResponse> => {
    return vulcanPost(vulcanPath(ROUTES.submit(workflowName)), session);
}

export const getData = (url: string) => {
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
