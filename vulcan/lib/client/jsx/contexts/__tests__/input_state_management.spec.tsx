import {integrateElement} from "../../test_utils/integration";
import {SessionStatusResponse, VulcanSession, WorkflowsResponse} from "../../api_types";
import {workflowsResponse} from "../../test_utils/fixtures/workflows-response";
import {act} from "react-test-renderer";
import {setInputs, setStatus, setWorkflow} from "../../actions/vulcan";
import {
  createStatusFixture,
  createStatusResponseFixture,
  createStepFixture, createStepStatusFixture,
  createWorkflowFixture,
  findWorkflowFromResponse
} from "../../test_utils/fixtures";
import {defaultContext} from "../vulcan_context";

describe('useInputStateManagement', () => {
  it('works', async () => {
    const {contextData, setInput, setData} = integrateElement(() => null, {
      providerOverrides: {
        getWorkflows: defaultContext.getWorkflows,
        pollStatus: defaultContext.pollStatus,
      }
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
      contextData.dispatch(setWorkflow(workflow));
      contextData.dispatch(setStatus(createStatusFixture(workflow)))
    });

    expect(contextData.state.inputs).toEqual({
      "query1/result": null,
    });

    await act(async function () {
      setInput('a', 10);
      setInput('b', 20);
    });

    expect(contextData.state.inputs).toEqual({
      "a": 10,
      "b": 20,
      "query1/result": null,
    });

    await act(async function () {
      setInput('query1/result', 30);
    });

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
    });
  })
})