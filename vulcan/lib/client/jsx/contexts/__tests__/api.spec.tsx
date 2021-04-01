import {act} from "react-test-renderer";
import * as React from "react";
import {stubUrl} from "etna-js/spec/helpers";
import {defaultVulcanSession} from "../../api_types";
import {workflowsResponse} from "../../test_utils/fixtures/workflows-response";
import {statusWithDownloads} from "../../test_utils/fixtures/status-with-downloads";
import {statusWithoutDownloads} from "../../test_utils/fixtures/status-without-downloads";
import {integrateElement} from "../../test_utils/integration";


describe('vulcanContext/useApiContextData', () => {
  it('works', async () => {
    const {contextData, reduxState, updateMatching} = integrateElement(() => null);

    await act(async function () {
      // The loading request fails, and results in an error being set.
      await updateMatching(() => reduxState.messages?.find((m: string) => m.indexOf('Our request was refused') !== -1));

      const expectedRequests = stubUrl({
          verb: 'get',
          url: 'https://vulcan.test/api/workflows',
          response: workflowsResponse,
        }).then(() => stubUrl({
        verb: 'get',
        url: 'https://download1',
        response: JSON.stringify({a: 1}),
      }).then(() => stubUrl({
        verb: 'get',
        url: 'https://download2',
        response(a, b, cb) {
          cb(null, [200, "1-23"]);
        }
      })).then(() => stubUrl({
        verb: 'post',
        url: 'https://vulcan.test/api/session//status',
        request: defaultVulcanSession,
        response: statusWithoutDownloads
      })).then(() => stubUrl({
        verb: 'post',
        url: 'https://vulcan.test/api/session/',
        request: defaultVulcanSession,
        response: statusWithDownloads
      })));

      await contextData.getWorkflows().then(r => expect(r).toEqual(workflowsResponse));
      await contextData.getData('https://download1').then(r => expect(r).toEqual({a: 1}));
      await contextData.getData('https://download2').then(r => expect(r).toEqual('1-23'));
      await contextData.pollStatus(defaultVulcanSession).then(r => expect(r).toEqual(statusWithoutDownloads));
      await contextData.postInputs(defaultVulcanSession).then(r => expect(r).toEqual(statusWithDownloads));

      await expectedRequests;

      expect(contextData.isLoading).toEqual(false);

      await Promise.all([
        updateMatching(() => contextData.isLoading),
        contextData.scheduleWork(Promise.resolve()),
      ])

      expect(contextData.isLoading).toEqual(false);

      await Promise.all([
        updateMatching(() => contextData.isLoading),
        await contextData.scheduleWork(Promise.reject(new Error('Oh snap!'))).catch(() => 0),
      ]);

      expect(reduxState.messages).toEqual(["Error: Oh snap!"]);
      expect(contextData.isLoading).toEqual(false);
    });

  })
});