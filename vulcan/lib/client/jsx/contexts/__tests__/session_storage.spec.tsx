import * as React from 'react';
import {
  awaitBefore,
  integrateElement,
  setupBefore
} from '../../test_utils/integration';
import {createFakeStorage} from '../../test_utils/mocks';
import {setSession, setSessionAndFigure} from '../../actions/vulcan_actions';
import {useWorkflowUtils} from '../../test_utils/workflow_utils';
import {useContext} from 'react';
import {VulcanContext} from '../vulcan_context';
import {
  createSessionFixture,
  createStoredSessionFixture
} from '../../test_utils/fixtures';

describe('useSessionStorage', () => {
  const integrated = setupBefore(integrateElement);
  const storage = setupBefore(() => {
    const storage = createFakeStorage();
    integrated.value.provideOverrides({
      storage
    });
    return storage;
  });
  const workflowHelpers = setupBefore(() =>
    integrated.value.runHook(() => useWorkflowUtils())
  );
  const contextData = setupBefore(() =>
    integrated.value.runHook(() => useContext(VulcanContext))
  );

  const beforeSessionUpdate = awaitBefore(async () => {});
  const nextSession = setupBefore(() =>
    createStoredSessionFixture('test', {project_name: 'test'})
  );
  const sessionUpdate = awaitBefore(async () => {
    await beforeSessionUpdate.ensure();
    await nextSession.ensure();
    contextData.value.dispatch(setSession(nextSession.value));
  });

  const sessionKeys = setupBefore(() => {
    const result: string[] = [];
    for (let i = 0; i < storage.value.length; ++i) {
      result.push(storage.value.key(i) || '');
    }
    return result;
  });

  const reloadedSession = awaitBefore(async () => {
    const {getLocalSession} = contextData.value;
    return await getLocalSession(
      workflowHelpers.value.workflow.name,
      nextSession.value.project_name,
      null
    );
  });

  it('does not save the session without a workflow being set', () => {
    expect(storage.value.length).toEqual(0);
  });

  describe('after having the workflow set', () => {
    const setupWorkflow = awaitBefore(async () => {
      workflowHelpers.value.setWorkflow('test');
    });

    beforeSessionUpdate.replace(async (factory) => {
      await factory();
      await setupWorkflow.ensure();
    });

    it('stores the session', async () => {
      expect(sessionKeys.value).toEqual([
        '/figure/new/',
        'test/figure/new/test'
      ]);
      expect(reloadedSession.value).toEqual(nextSession.value);
    });

    describe('after an secondary session update', () => {
      const previousSession = setupBefore(() =>
        createSessionFixture('test', {project_name: 'test'})
      );
      nextSession.replace(async (factory) => {
        const original = await factory();
        return {...original, key: 'abcdef'};
      });

      beforeSessionUpdate.replace(async (factory) => {
        await nextSession.ensure();
        await factory();
        await previousSession.ensure();
        contextData.value.dispatch(setSession(previousSession.value));
      });

      it('updates that value', async () => {
        expect(reloadedSession.value).toEqual({
          figure_id: null,
          inputs: {},
          key: 'abcdef',
          project_name: 'test',
          workflow_name: 'test'
        });
      });

      describe('with a different project_name', () => {
        beforeSessionUpdate.replace(async (factory) => {
          await factory();
          workflowHelpers.value.setWorkflow('test', {}, ['test2'], 'test2');
        });

        nextSession.replace(async (factory) => {
          const original = await factory();
          return {...original, project_name: 'test2'};
        });

        it('stores that session in a separate key', () => {
          expect(reloadedSession.value).toEqual({
            figure_id: null,
            inputs: {},
            key: 'abcdef',
            project_name: 'test2',
            workflow_name: 'test'
          });
        });
      });
    });

    describe('with figure information', () => {
      const previousStoredSession = setupBefore(() =>
        createStoredSessionFixture(
          'test',
          {
            key: 'abcdef',
            project_name: 'test3',
            inputs: {a: '123', c: 'def'}
          },
          {
            figure_id: 1,
            inputs: {a: '123', c: 'def'}
          }
        )
      );
      const userSession = setupBefore(() =>
        createSessionFixture('test', {
          project_name: 'test3',
          key: 'session-key',
          inputs: {a: '321', b: 'abdef'}
        })
      );
      reloadedSession.replace(async (factory) => {
        const {getLocalSession} = contextData.value;
        return await getLocalSession(
          workflowHelpers.value.workflow.name,
          userSession.value.project_name,
          1
        );
      });

      beforeSessionUpdate.replace(async (factory) => {
        await userSession.ensure();
        await factory();
        await previousStoredSession.ensure();
        contextData.value.dispatch(
          setSessionAndFigure(previousStoredSession.value)
        );
        contextData.value.dispatch(setSession(userSession.value));
      });

      it('stores updated session inputs', () => {
        expect(reloadedSession.value).toEqual({
          figure_id: 1,
          inputs: {
            a: '321',
            b: 'abdef'
          },
          key: 'session-key',
          project_name: 'test3',
          workflow_name: 'test'
        });
      });
    });
  });

  it('works', async () => {
    // const {dispatch, stateRef} = contextData.value;
    // const storage = createFakeStorage();
    // await dispatch(setSession(statusWithDownloads['session']));
    // //
    // // Does not save without workflow being set
    // expect(storage.length).toEqual(0);
    // await dispatch(setWorkflow({...defaultWorkflow, name: 'some name', projects: ['a project']}, 'a project'));
    //
    // expect(storage.length).toEqual(1); // Saved it!
    // expect(contextData.state.session.workflow_name).toEqual(contextData.state.workflow?.name)
    // expect(await contextData.getLocalSession(contextData.state.workflow as any, 'a project')).toEqual(contextData.state.session);
    // await dispatch(setSession({
    //   ...statusWithDownloads['session'],
    //   key: 'anewkey',
    //   project_name: 'a project',
    //   workflow_name: 'some name',
    // }));
    // expect(contextData.state.session.key).toEqual('anewkey');
    // expect(contextData.state.session.project_name).toEqual('a project');
    // expect(storage.length).toEqual(1); // Updates it!
    // expect(contextData.state.workflow).not.toBeFalsy();
    //
    // expect(await contextData.getLocalSession(contextData.state.workflow as any, 'a project')).toEqual(contextData.state.session);
    //
    // await dispatch(setWorkflow({...defaultWorkflow, name: 'new name', projects: ['test', 'test2'] }, 'test'));
    // await dispatch(setSession({
    //   ...statusWithDownloads['session'],
    //   workflow_name: 'new name',
    // }));
    //
    // expect(storage.length).toEqual(2); // Adds the new work flow session key
    //
    // // Switch projects while retaining the same workflow.
    // await dispatch(setWorkflow(contextData.stateRef.current.workflow as Workflow, 'test2'));
    // expect(storage.length).toEqual(3); // Adds the new work flow session key
  });
});
