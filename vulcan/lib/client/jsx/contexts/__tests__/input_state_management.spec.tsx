import {integrateElement} from "../../test_utils/integration";
import {act} from "react-test-renderer";
import {patchInputs, setStatus, setWorkflow} from "../../actions/vulcan";
import {
  createStatusFixture,
  createStepFixture,
  createWorkflowFixture,
} from "../../test_utils/fixtures";
import {defaultContext} from "../vulcan_context";

describe('useInputStateManagement', () => {
  it('works', async () => {
    const overrides = {
      getWorkflows: defaultContext.getWorkflows,
      pollStatus: defaultContext.pollStatus,
      statusIsFresh: true,
    };

    const {contextData, dispatch, setData, replaceOverrides} = integrateElement(() => null, {
      providerOverrides: overrides
    });

    const workflow = createWorkflowFixture({
      inputs: {
        a: {type: 'int'},
        b: {type: 'int'},
      },

      dependencies_of_outputs: {
        'a': ['process/result', 'query2/result'],
        'b': [],
        'query1/result': ['process/result', 'query2/result'],
        'process/result': ['query2/result'],
        'query2/result': [],
      },

      steps: [
        [
          createStepFixture({name: 'query1', run: 'ui-queries/thing', out: ['result']}),
          createStepFixture({
            name: 'process',
            run: 'scripts/process',
            in: [{source: 'a', id: 'a'}, {source: 'query1/result', id: 'b'}],
            out: ['result']
          }),
          createStepFixture({
            name: 'query2',
            run: 'ui-queries/thing',
            in: [{source: 'process/result', id: 'a'}],
            out: ['result']
          }),
        ]
      ]
    })

    await act(async function () {
      contextData.dispatch(setWorkflow(workflow, 'test'));
      contextData.dispatch(setStatus(createStatusFixture(workflow)))
    });

    expect(contextData.state.inputs).toEqual({});

    await dispatch(patchInputs({ a: 10, b: 20 }));

    expect(contextData.state.inputs).toEqual({
      "a": 10,
      "b": 20,
      "query1/result": null,
    });

    await dispatch(patchInputs({ 'query1/result': 30 }));

    expect(contextData.state.inputs).toEqual({
      "a": 10,
      "b": 20,
      "query1/result": 30,
    });

    await act(async function () {
      setData('process/result', 30);
    });

    expect(contextData.state.inputs).toEqual({
      "a": 10,
      "b": 20,
      "query1/result": 30,
      "query2/result": null,
    });

    expect(contextData.state.session.inputs).toEqual({
      "a": 10,
      "b": 20,
      "query1/result": 30,
    });

    await dispatch(patchInputs({ 'query1/result': null }));

    // query2 is dropped as a result of a dropping the download of process as a result of dropping query1.
    expect(contextData.state.inputs).toEqual({
      "a": 10,
      "b": 20,
      "query1/result": null,
    });

    expect(contextData.state.session.inputs).toEqual({
      "a": 10,
      "b": 20,
    });

    // Allow invalid states while status is till being fetched.
    replaceOverrides({...overrides, statusIsFresh: false});

    await dispatch(patchInputs({ 'query2/result': 33 }));
    expect(contextData.state.session.inputs).toEqual({
      "a": 10,
      "b": 20,
      "query2/result": 33,
    });
  })
})