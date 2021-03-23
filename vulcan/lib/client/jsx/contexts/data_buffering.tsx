import {Dispatch, useEffect, useState} from "react";
import {releaseDownloadedData, setDownloadedData, VulcanAction} from "../actions/vulcan";
import {defaultApiHelpers} from "./api";
import {getData} from "../api/vulcan_api";
import {VulcanState} from "../reducers/vulcan_reducer";
import {shouldDownload, stepOfStatus} from "../selectors/workflow_selectors";

export function useDataBuffering(
    state: VulcanState,
    dispatch: Dispatch<VulcanAction>,
    scheduleWork: typeof defaultApiHelpers.scheduleWork
) {
    const [activeDownloads, setActiveDownloads] = useState({} as {[k: string]: boolean});

    useEffect(() => {
        // Downloads we found in a status.
        const foundDownloads: {[k: string]: boolean} = {};
        const {workflow} = state;

        if (workflow != null) {
            state.status[0].forEach(stepStatus => {
                if (!stepStatus.downloads) return;

                Object.values(stepStatus.downloads).forEach(url => {
                    if (!shouldDownload(url, workflow, stepOfStatus(stepStatus, workflow), state.status)) return;

                    foundDownloads[url] = true;
                    if (activeDownloads[url] || url in state.data) return;

                    // Mark that we are trying to fetch the data already, to prevent duplicate downloads.
                    setActiveDownloads({...activeDownloads, [url]: true});

                    scheduleWork(getData(url)).then(data => {
                        dispatch(setDownloadedData(url, data));
                        setActiveDownloads({...activeDownloads, [url]: false});
                    });
                })
            });
        }

        Object.keys(state.data).forEach(url => {
            if (url in foundDownloads) return;
            dispatch(releaseDownloadedData(url));
        });
    }, [state.data, state.status, state.workflow]);
}
