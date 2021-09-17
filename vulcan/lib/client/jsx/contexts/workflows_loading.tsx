import {Dispatch, useEffect} from "react";
import {setWorkflows, VulcanAction} from "../actions/vulcan_actions";
import {defaultApiHelpers} from "./api";
import {Cancellable} from "etna-js/utils/cancellable";

export function useWorkflowsLoading(
    key: string,
    dispatch: Dispatch<VulcanAction>,
    getWorkflows: typeof defaultApiHelpers.getWorkflows,
    showErrors: typeof defaultApiHelpers.showErrors) {

    useEffect(() => {
        const cancellable = new Cancellable();

        showErrors(cancellable.race(getWorkflows())
            .then(({result, cancelled}) => {
                if (result && !cancelled) dispatch(setWorkflows(result.workflows));
            }));

        return () => cancellable.cancel();
    }, [key, dispatch, showErrors, getWorkflows]);
}