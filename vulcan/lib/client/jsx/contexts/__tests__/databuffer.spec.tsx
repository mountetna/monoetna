import {awaitBefore, integrateElement, setupBefore} from "../../test_utils/integration";
import {VulcanContext} from "../vulcan_context";
import {useWorkflowUtils} from "../../test_utils/workflow_utils";
import {useContext} from "react";
import {act} from "react-test-renderer";

describe('useDataBuffering', () => {
  beforeEach(() => {
    jest.useFakeTimers();
  })

  afterEach(() => {
    jest.useRealTimers();
  })

  const integrated = setupBefore(integrateElement);
  const getDataMock = setupBefore(() => integrated.value.blockingAsyncMock('getData'));
  const workflowHelpers = setupBefore(() => integrated.value.runHook(() => useWorkflowUtils()));
  const contextData = setupBefore(() => integrated.value.runHook(() => useContext(VulcanContext)));

  const setupWorkflow = awaitBefore(async () => {
    workflowHelpers.value.setWorkflow('test');
  })

  const aStep = awaitBefore(async () => {
    return workflowHelpers.value.addStep('astep');
  })

  const aStepStatus = awaitBefore(async () => {
    return workflowHelpers.value.setStatus('astep', {
      status: 'complete',
      downloads: { 'a': 'http://astep-a', }
    });
  })

  const dependent = awaitBefore(async () => {
    return workflowHelpers.value.addStep('uiQueryStep', {
      run: 'ui-queries/something',
      in: [{ id: 'cool', source: aStep.value.name + '/a' }]
    })
  })

  it('downloads data from that url', async () => {
    const {stateRef} = contextData.value;
    await getDataMock.value.respond(async (url) => {
      expect(url).toEqual("http://astep-a");
      return "somedata";
    })

    expect(stateRef.current.data).toEqual({
      "http://astep-a": "somedata",
    });
    expect(getDataMock.value.pendingCount()).toEqual(0);
  })

  describe('with multiple downloads', () => {
    aStepStatus.replace(async factory => {
      const updated = await factory();
      updated.downloads = {...updated.downloads, b: 'http://astep-b'};
      workflowHelpers.value.setStatus('astep', updated);
      return updated;
    })

    const bStep = awaitBefore(async () => {
      return workflowHelpers.value.addStep('bstep');
    });

    const bStepStatus = awaitBefore(async () => {
      return workflowHelpers.value.setStatus('bstep', {
        status: 'complete',
        downloads: { 'a': 'http://bstep-a', }
      });
    })

    dependent.replace(async () => {
      await bStep.ensure();
      await aStep.ensure();
      return workflowHelpers.value.addStep('uiQueryStep', {
        run: 'ui-queries/something',
        in: [
          { id: 'cool', source: aStep.value.name + '/a' },
          { id: 'cool', source: aStep.value.name + '/b' },
          { id: 'cool', source: bStep.value.name + '/a' },
        ]
      });
    })

    const afterTheFirstDownload = awaitBefore(async () => {
      await getDataMock.value.respond(async (url) => {
        expect(url).toEqual("http://astep-a");
        return "thing1";
      });
    });

    describe('when there is a failure', () => {
      afterTheFirstDownload.replace(async () => {
        await getDataMock.value.reject((url) => {
          throw new Error('Mock Request Failure');
        });
      })

      it('retries, then continues downloading', async () => {
        const {stateRef} = contextData.value;
        expect(getDataMock.value.pendingCount()).toEqual(0);

        await act(async () => { jest.advanceTimersByTime(5000) });
        expect(getDataMock.value.pendingCount()).toEqual(1);

        await getDataMock.value.respond(async (url) => {
          expect(url).toEqual("http://astep-a");
          return 1;
        });

        await getDataMock.value.respond(async (url) => {
          expect(url).toEqual("http://astep-b");
          return 1;
        });

        expect(getDataMock.value.pendingCount()).toEqual(1);
      })
    })

    it('begins and continues to download the remaining', async () => {
      const {stateRef} = contextData.value;
      expect(getDataMock.value.pendingCount()).toEqual(1);

      await getDataMock.value.respond(async (url) => {
        expect(url).toEqual("http://astep-b");
        return 2;
      });

      expect(getDataMock.value.pendingCount()).toEqual(1);
      await getDataMock.value.respond(async (url) => {
        expect(url).toEqual("http://bstep-a");
        return 3;
      });

      expect(getDataMock.value.pendingCount()).toEqual(0);
      expect(stateRef.current.data).toEqual({
        "http://astep-a": "thing1",
        "http://astep-b": 2,
        "http://bstep-a": 3,
      });
    })
  })

  describe('when the dependency does not reference the download directly', () => {
    dependent.replace(async (factory) => {
      const original = await factory();
      original.in[0].source = `${aStep.value.name}/b`;
      return original;
    })

    it('downloads nothing', async () => {
      expect(getDataMock.value.pendingCount()).toEqual(0);
    })
  })

  describe('when the dependency is a script', () => {
    dependent.replace(async (factory) => {
      const original = await factory();
      original.run = 'script/do-a-thing'
      return original;
    })

    it('downloads nothing', async () => {
      expect(getDataMock.value.pendingCount()).toEqual(0);
    })
  })

  describe('when the dependency is a link', () => {
    dependent.replace(async (factory) => {
      const original = await factory();
      original.run = 'ui-outputs/link'
      return original;
    })

    it('downloads nothing', async () => {
      expect(getDataMock.value.pendingCount()).toEqual(0);
    })
  })
});