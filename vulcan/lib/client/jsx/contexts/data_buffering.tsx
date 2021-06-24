import {Dispatch, useEffect, useMemo, useState} from "react";
import {releaseDownloadedData, setDownloadedData, VulcanAction} from "../actions/vulcan_actions";
import {defaultApiHelpers} from "./api";
import {VulcanState} from "../reducers/vulcan_reducer";
import {shouldDownload, stepOfStatus} from "../selectors/workflow_selectors";

export const defaultDataBufferingHelpers = {
  urlsToBuffer: {} as { [k: string]: boolean },
  missingDownloads: [] as string[],
  activeDownload: null as string | null,
}

export function useDataBuffering(
    state: VulcanState,
    dispatch: Dispatch<VulcanAction>,
    scheduleWork: typeof defaultApiHelpers.scheduleWork,
    getData: typeof defaultApiHelpers.getData,
): typeof defaultDataBufferingHelpers {
  const [activeDownload, setActiveDownload] = useState(null as string | null);
  const {status, workflow, data} = state;

  const urlsToBuffer = useMemo(() => {
    const result: { [k: string]: boolean } = {};
    if (workflow != null) {
      status[0].forEach(stepStatus => {
        if (!stepStatus.downloads) return;
        Object.values(stepStatus.downloads).forEach(url => {
          if (!shouldDownload(url, workflow, stepOfStatus(stepStatus, workflow), status)) return;
          result[url] = true;
        });
      });
    }

    return result;
  }, [status, workflow])

  const missingDownloads = useMemo(() => {
    return Object.keys(urlsToBuffer).filter(url => {
      return !(activeDownload == url || url in data);
    });
  }, [urlsToBuffer, data, activeDownload]);

  // Queue up and send out download requests
  useEffect(() => {
    const nextDownload = missingDownloads[0];
    if (!nextDownload) return;
    if (activeDownload != null) return;

    setActiveDownload(nextDownload);

    scheduleWork(getData(nextDownload).then(data => {
      dispatch(setDownloadedData(nextDownload, data));
    })).catch(e => { console.error(e); }).finally(() => {
      setActiveDownload(null);
    });
  }, [missingDownloads, activeDownload, scheduleWork, getData, dispatch]);

  // Clear out unused downloads, one key at a time.
  useEffect(() => {
    const toRelease = Object.keys(data).find(url => !(url in urlsToBuffer));
    if (toRelease != null) {
      dispatch(releaseDownloadedData(toRelease));
    }
  }, [urlsToBuffer, data, dispatch])

  return {
    urlsToBuffer,
    missingDownloads,
    activeDownload,
  };
}
