import {createFakeStorage} from "../../test_utils/mocks";
import {integrateElement} from "../../test_utils/integration";
import {act} from "react-test-renderer";
import {setStatus, setWorkflow} from "../../actions/vulcan";
import {createStatusFixture, createStepStatusFixture, findWorkflowFromResponse} from "../../test_utils/fixtures";
import {workflowsResponse} from "../../test_utils/fixtures/workflows-response";
import {defaultContext} from "../vulcan_context";
import {delay} from "etna-js/spec/helpers";
import {defaultWorkflow} from "../../api_types";

describe('useDataBuffering', () => {
  it('eventually downloads all the things, and clears when necessary', async () => {
    const storage = createFakeStorage();
    const q: [(v: any) => void, (e: any) => void][] = [];

    async function answerDownload(v: any) {
      await act(async function() {
        const next = q.shift();
        if (next) {
          next[0](v);
          await updateMatching(() => contextData.activeDownload != null);
        }
      })
    }

    async function failDownload(e: any) {
      await act(async function() {
        const next = q.shift();
        if (next) {
          next[1](e);
          await updateMatching(() => contextData.activeDownload != null);
        }
      });
    }

    const {contextData, updateMatching, dispatch} = integrateElement(() => null,
        {
          providerOverrides: {
            getWorkflows: defaultContext.getWorkflows,
            storage, getData(url: string): Promise<any> {
              return new Promise<any>((resolve, reject) => {
                q.push([resolve, reject]);
              })
            }
          }
        });

    let workflow = findWorkflowFromResponse(workflowsResponse, "test_concurrent_workflow.cwl");

    await dispatch(setWorkflow(workflow));
    await dispatch(setStatus(createStatusFixture(workflow,
        createStepStatusFixture({name: 'firstAdd', downloads: {'sum': 'https://firstAddSum'}}),
        createStepStatusFixture({
          name: 'twiceOtherPickANumber',
          downloads: {'sum': 'https://twiceOtherPickANumberSum'}
        }),
        createStepStatusFixture({
          name: 'twiceOtherPickANumber2',
          downloads: {'sum': 'https://twiceOtherPickANumber2Sum'}
        }),
        createStepStatusFixture({
          name: 'thriceOtherPickANumber',
          downloads: {'sum': 'https://thriceOtherPickANumberSum'}
        }),
        createStepStatusFixture({name: 'ficeOtherPickANumber', downloads: {'sum': 'https://ficeOtherPickANumberSum'}}),
        createStepStatusFixture({name: 'finalStep', downloads: {'sum': 'https://finalStepSum'}}),
    )));

    const originalUrlsToBuffer = {...contextData.urlsToBuffer};
    expect(contextData.urlsToBuffer).toEqual({
      "https://finalStepSum": true,
      "https://firstAddSum": true,
    });

    function checkState(expectedMissing: number, shouldHaveActiveDownload = true) {
      expect(contextData.urlsToBuffer).toEqual(originalUrlsToBuffer);

      if (shouldHaveActiveDownload) expect(Object.keys(contextData.urlsToBuffer)).toContain(contextData.activeDownload);
      expect(contextData.missingDownloads.length).toEqual(expectedMissing);
      contextData.missingDownloads.forEach(k => expect(Object.keys(contextData.urlsToBuffer)).toContain(k));
      contextData.missingDownloads.forEach(k => expect(contextData.missingDownloads.filter(v => v === k).length).toEqual(1));
      expect(contextData.missingDownloads).not.toContain(contextData.activeDownload);
      if (shouldHaveActiveDownload) expect(Object.keys(contextData.urlsToBuffer)).toContain(contextData.activeDownload);
    }

    checkState(1);
    await failDownload(new Error("Oh no!"));
    checkState(1);
    await answerDownload(123);
    checkState(0);
    await answerDownload(700);
    checkState(0, false);

    expect(Object.keys(contextData.state.data).sort()).toEqual(Object.keys(contextData.urlsToBuffer).sort());
    expect(Object.values(contextData.state.data).sort()).toEqual([123, 700]);

    expect(contextData.urlsToBuffer).toEqual(originalUrlsToBuffer);

    await delay(1000);
    // No extra, spare requests sitting around.
    expect(q.length).toEqual(0);

    // Clears data when it is not needed
    await dispatch(setWorkflow(defaultWorkflow));
    expect(contextData.urlsToBuffer).toEqual({});
    expect(contextData.state.data).toEqual({});
    expect(contextData.activeDownload).toEqual(null);
    expect(q.length).toEqual(0);
  });
});