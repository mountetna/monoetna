import React, {useCallback, useContext, useEffect, useMemo} from 'react';
import Icon from 'etna-js/components/icon';
import Link from 'etna-js/components/link';

import {VulcanContext} from '../../../contexts/vulcan_context';
import InputFeed from './input_feed';
import OutputFeed from './output_feed';
import {statusOfStep, uiOutputOfStep, workflowName} from "../../../selectors/workflow_selectors";
import {commitInputs} from "../../../actions/vulcan";

export default function SessionManager() {
  const context = useContext(VulcanContext);
  const {state, dispatch} = context;

  const workflow = state.workflow;
  if (!workflow) return null;
  const name = workflowName(workflow)
  if (!name) return null;

  const {steps} = workflow;
  const {status} = state;

  // We are done once every step either has a download or that step is a uiOutput.
  const complete = useMemo(() =>
          steps[0].every(step => uiOutputOfStep(step) || statusOfStep(step, status)?.downloads),
      [steps, status])

  const run = useCallback(() => {
    dispatch(commitInputs());
  }, [dispatch, commitInputs])

  return (
    <div className='session-manager'>
      <div className='session-header'>
        <span className='session-workflow-name'>{name}</span>
        {workflow.vignette && (
          <Link link={ROUTES.workflow_vignette(name)}>
            <Icon className='vignette' icon='book' />
          </Link>
        )}
        <Icon
          className='run'
          disabled={!complete}
          title='Run workflow'
          onClick={run}
          icon='play'
        />
      </div>
      <div className='session-feed-container'>
        <InputFeed/>
        <OutputFeed/>
      </div>
    </div>
  );
}
