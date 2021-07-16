import {act} from "react-test-renderer";
import * as React from "react";
import {integrateElement} from "../../test_utils/integration";
import {createFakeStorage} from "../../test_utils/mocks";
import {setSession, setWorkflow} from "../../actions/vulcan_actions";
import {statusWithDownloads} from "../../test_utils/fixtures/status-with-downloads";
import {defaultWorkflow, Workflow} from "../../api_types";

describe('useSessionStorage', () => {
  it('works', async () => {
    const storage = createFakeStorage();

    const {contextData, dispatch} = integrateElement(() => null, {
      providerOverrides: {
        storage
      }
    });

    await dispatch(setSession(statusWithDownloads['session']));

    // Does not save without workflow being set
    expect(storage.length).toEqual(0);
    await dispatch(setWorkflow({...defaultWorkflow, name: 'some name', projects: ['a project']}, 'a project'));

    expect(storage.length).toEqual(1); // Saved it!
    expect(contextData.state.session.workflow_name).toEqual(contextData.state.workflow?.name)
    expect(await contextData.getLocalSession(contextData.state.workflow as any, 'a project')).toEqual(contextData.state.session);

    await dispatch(setSession({
      ...statusWithDownloads['session'],
      key: 'anewkey',
      project_name: 'a project',
      workflow_name: 'some name',
    }));
    expect(contextData.state.session.key).toEqual('anewkey');
    expect(contextData.state.session.project_name).toEqual('a project');
    expect(storage.length).toEqual(1); // Updates it!
    expect(contextData.state.workflow).not.toBeFalsy();

    expect(await contextData.getLocalSession(contextData.state.workflow as any, 'a project')).toEqual(contextData.state.session);

    await dispatch(setWorkflow({...defaultWorkflow, name: 'new name', projects: ['test', 'test2'] }, 'test'));
    await dispatch(setSession({
      ...statusWithDownloads['session'],
      workflow_name: 'new name',
    }));

    expect(storage.length).toEqual(2); // Adds the new work flow session key

    // Switch projects while retaining the same workflow.
    await dispatch(setWorkflow(contextData.stateRef.current.workflow as Workflow, 'test2'));
    expect(storage.length).toEqual(3); // Adds the new work flow session key
  });
});
