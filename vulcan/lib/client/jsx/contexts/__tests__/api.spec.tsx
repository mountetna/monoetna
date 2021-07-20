import {act} from "react-test-renderer";
import * as React from "react";
import {stubUrl} from "etna-js/spec/helpers";
import {defaultVulcanSession} from "../../api_types";
import {workflowsResponse} from "../../test_utils/fixtures/workflows-response";
import {statusWithDownloads} from "../../test_utils/fixtures/status-with-downloads";
import {statusWithoutDownloads} from "../../test_utils/fixtures/status-without-downloads";
import {integrateElement, setupBefore} from "../../test_utils/integration";
import {useContext} from "react";
import {VulcanContext} from "../vulcan_context";

describe('contexts/api', () => {
  const integrated = setupBefore(() => integrateElement(null, {
    getWorkflows: undefined,
    getData: undefined,
    pollStatus: undefined,
    postInputs: undefined,
    showErrors: undefined,
  }));
  const contextData = setupBefore(() => integrated.value.runHook(() => useContext(VulcanContext)));

  it('works', async () => {
    const { getWorkflows, getData, pollStatus, postInputs, showErrors } = contextData.value;

    const session = {...defaultVulcanSession, workflow_name: 'abc', project_name: 'test'};

    await act(async function () {
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
        url: 'https://vulcan.test/api/test/session/abc/status',
        request: session,
        response: statusWithoutDownloads
      })).then(() => stubUrl({
        verb: 'post',
        url: 'https://vulcan.test/api/test/session/abc',
        request: session,
        response: statusWithDownloads
      })));

      await getWorkflows().then(r => expect(r).toEqual(workflowsResponse));
      await getData('https://download1').then(r => expect(r).toEqual({a: 1}));
      await getData('https://download2').then(r => expect(r).toEqual('1-23'));
      await pollStatus(session).then(r => expect(r).toEqual(statusWithoutDownloads));
      await postInputs(session).then(r => expect(r).toEqual(statusWithDownloads));

      await expectedRequests;

      await showErrors(Promise.reject(new Error('Oh snap!'))).catch(() => 0);
    });

  })
});