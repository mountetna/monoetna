import {Dispatch, useEffect} from "react";
import {setWorkflows, VulcanAction} from "../actions/vulcan";
import {defaultApiHelpers} from "./api";
import {Cancellable} from "etna-js/utils/cancellable";

export function useWorkflowsLoading(
    key: string,
    dispatch: Dispatch<VulcanAction>,
    getWorkflows: typeof defaultApiHelpers.getWorkflows,
    scheduleWork: typeof defaultApiHelpers.scheduleWork) {

    useEffect(() => {
        const cancellable = new Cancellable();

        scheduleWork(cancellable.race(getWorkflows())
            .then(({result, cancelled}) => {
                if (result && !cancelled) setWorkflows(result.workflows);
            }));

        return () => cancellable.cancel();
    }, [key]);
}