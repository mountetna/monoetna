import {Dispatch} from "react";
import {setDownloadedData, VulcanAction} from "../actions/vulcan_actions";
import {defaultApiHelpers} from "./api";
import {VulcanState} from "../reducers/vulcan_reducer";
import {shouldDownloadStep} from "../selectors/workflow_selectors";
import {useAsync} from "etna-js/utils/cancellable_helpers";
import {runAttempts} from "etna-js/utils/retryable";

export function useDataBuffering(
    state: VulcanState,
    dispatch: Dispatch<VulcanAction>,
    showErrors: typeof defaultApiHelpers.showErrors,
    getData: typeof defaultApiHelpers.getData,
) {
  const {status, workflow, data} = state;

  useAsync(function* () {
    if (!workflow) return;

    // Find the next status with a download that is an input to a ui-output or ui-query,
    // initiate a download, and let it run.
    for (let stepStatus of status[0]) {
      if (!stepStatus.downloads) continue;
      for (let download of Object.entries(stepStatus.downloads)) {
        const [outputName, url] = download;
        if (url in data) continue;
        if (!shouldDownloadStep(stepStatus.name, workflow, outputName)) continue;

        const downloaded = yield* runAttempts(() => getData(url));

        dispatch(setDownloadedData(url, downloaded));
        return;
      }
    }
  }, [status, workflow, data]);
}
