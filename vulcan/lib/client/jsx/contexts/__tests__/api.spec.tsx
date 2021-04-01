import {act, create} from "react-test-renderer";
import {defaultContext, VulcanContext, VulcanContextData, VulcanProvider} from "../vulcan_context";
import {useContext} from "react";
import {Provider} from "react-redux";
import * as React from "react";
import {VulcanStore} from "../../vulcan_store";
import {useReduxState} from "etna-js/hooks/useReduxState";
import {joinedDeferredPromises, stubUrl} from "etna-js/spec/helpers";
import {defaultVulcanSession, defaultWorkflowsResponse} from "../../api_types";
import {workflowsResponse} from "../../test_utils/fixtures/workflows-response";
import {statusWithDownloads} from "../../test_utils/fixtures/status-with-downloads";
import {statusWithoutDownloads} from "../../test_utils/fixtures/status-without-downloads";


describe('vulcanContext/useApiHelpers', () => {
  it('works', async () => {
    let helpers: VulcanContextData = defaultContext;
    let reduxState: any;
    let waiters: Function[] = [];

    function updateMatching(pred: () => boolean): Promise<void> {
      return new Promise((resolve, reject) => {
        let timedOut = false;

        const t = setTimeout(function () {
          timedOut = true;
          reject('timeout');
        }, 3000);

        waiters.push(function () {
          if (!pred()) return false;
          if (timedOut) return true;

          clearTimeout(t);
          resolve();
          return true;
        });
      })
    }

    function TestComponent() {
      Object.assign(helpers, useContext(VulcanContext));
      reduxState = useReduxState();
      waiters.forEach(f => f());

      return null;
    }

    const node = create(
        <Provider store={VulcanStore()}>
          <VulcanProvider>
            <TestComponent/>
          </VulcanProvider>
        </Provider>
    )

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

      await helpers.getWorkflows().then(r => expect(r).toEqual(workflowsResponse));
      await helpers.getData('https://download1').then(r => expect(r).toEqual({a: 1}));
      await helpers.getData('https://download2').then(r => expect(r).toEqual('1-23'));
      await helpers.pollStatus(defaultVulcanSession).then(r => expect(r).toEqual(statusWithoutDownloads));
      await helpers.postInputs(defaultVulcanSession).then(r => expect(r).toEqual(statusWithDownloads));

      await expectedRequests;

      expect(helpers.isLoading).toEqual(false);

      await Promise.all([
        updateMatching(() => helpers.isLoading),
        helpers.scheduleWork(Promise.resolve()),
      ])

      expect(helpers.isLoading).toEqual(false);

      await Promise.all([
        updateMatching(() => helpers.isLoading),
        await helpers.scheduleWork(Promise.reject(new Error('Oh snap!'))).catch(() => 0),
      ]);

      expect(reduxState.messages).toEqual(["Error: Oh snap!"]);
      expect(helpers.isLoading).toEqual(false);
    });

  })
});