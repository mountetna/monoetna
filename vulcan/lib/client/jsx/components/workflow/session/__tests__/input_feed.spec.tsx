import React from 'react';
import {defaultContext, VulcanProvider} from '../../../../contexts/vulcan_context';
import renderer from 'react-test-renderer';
import InputFeed from '../input_feed';
import {stateFromActions} from "../../../../test_utils/state";
import {
    createStatusFixture,
    createStepFixture,
    createStepStatusFixture,
    createWorkflowFixture,
} from "../../../../test_utils/fixtures";
import {setStatus, setWorkflow, setWorkflows} from "../../../../actions/vulcan";

describe('InputFeed', () => {
    it('renders complete UI steps and error steps', () => {
        const workflow = createWorkflowFixture({
            inputs: {}, steps: [
                [
                    createStepFixture({name: 'zero'}),
                    createStepFixture({name: 'first', out: ['output']}),
                    createStepFixture({
                        name: 'second',
                        run: 'ui-queries/ask-the-user.cwl',
                        in: [{id: 'a', source: 'first/output'}],
                        out: ['response']
                    }),
                    createStepFixture({
                        name: 'third',
                        run: 'ui-queries/ask-again.cwl',
                        in: [{id: 'a', source: 'first/output'}],
                        out: ['response']
                    }),
                    createStepFixture({
                        name: 'fourth',
                        run: 'ui-outputs/show-the-user.cwl',
                    }),
                ]
            ]
        });

        const {state} = stateFromActions([
            setWorkflows([workflow]),
            setWorkflow(workflow),
            setStatus(createStatusFixture(workflow,
                createStepStatusFixture({name: 'zero', status:'error', error: 'Ooops!'}),
                createStepStatusFixture({name: 'first', status: 'complete', downloads: {output: 'https://download1'}}),
            )),
        ]);

        const component = renderer.create(
            <VulcanProvider state={state} useActionInvoker={defaultContext.useActionInvoker}>
                <InputFeed/>
            </VulcanProvider>
        );

        let instance = component.root;

        expect(instance.findAllByProps({className: 'step-error'}).length).toEqual(
            1
        );

        expect(
            instance.findAllByProps({className: 'step-user-input step'}).length
        ).toEqual(2);

        expect(component.toJSON()).toMatchSnapshot();
    });
});
