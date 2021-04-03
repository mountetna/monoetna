import React, {useCallback, useContext, useEffect, useMemo} from 'react';
import Icon from 'etna-js/components/icon';
import Link from 'etna-js/components/link';

import {VulcanContext} from '../../../contexts/vulcan_context';
import InputFeed from './input_feed';
import OutputFeed from './output_feed';
import {
    allWorkflowPrimaryInputSources,
    statusOfStep,
    uiOutputOfStep,
    workflowName
} from "../../../selectors/workflow_selectors";

export default function SessionManager() {
    const context = useContext(VulcanContext);
    const {state, dispatch, requestPoll} = context;

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

    const idle = useMemo(() =>
            steps[0].every(step => uiOutputOfStep(step) || statusOfStep(step, status)?.status !== 'running'),
        [steps, status]);

    const primaryInputsReady = useMemo(
        () => allWorkflowPrimaryInputSources(workflow).every(source => source in state.inputs),
        [state.inputs, workflow])

    const running = !idle;

    const run = useCallback(() => {
        console.log("hello?")
        requestPoll(true);
    }, [requestPoll])

    return (
        <div className='session-manager'>
            <div className='session-header'>
                <span className='session-workflow-name'>
                    {workflow.description || name}
                </span>
                {workflow.vignette && (
                    <div className='header-btn'>
                        <Link link={ROUTES.workflow_vignette(name)}>
                            Vignette
                            <Icon className='vignette' icon='book' />
                        </Link>
                    </div>
                )}
                <button
                    disabled={complete || running || !primaryInputsReady}
                    onClick={run}
                    className='run-workflow-btn header-btn'
                >
                    Run
                    <Icon
                        className='run'
                        disabled={complete || running || !primaryInputsReady}
                        title='Run workflow'
                        icon='play'
                    />
                </button>
            </div>
            <div className='session-feed-container'>
                <InputFeed/>
                <OutputFeed/>
            </div>
        </div>
    );
}
