import React from 'react';
import { rest } from 'msw';
import { setupServer } from 'msw/node';
import { fireEvent, waitFor, waitForElementToBeRemoved, screen } from '@testing-library/react';
import _ from 'lodash';

import NamesCreate from '../names-create';
import { renderWithProviders } from '../../../utils/test';
import * as exampleRulesResponse from './static/example-rules-response.json';
import { MagmaRulesResponse } from '../../../utils/rules';
import { MagmaBulkGenerateName, MagmaListName } from '../../../utils/names';
import * as exportModule from '../../../utils/export';
import { selectCompleteCreateNamesCreationPayloads } from '../../../selectors/names';


global.ResizeObserver = jest.fn().mockImplementation(() => ({
    observe: jest.fn(),
    unobserve: jest.fn(),
    disconnect: jest.fn(),
}));

const magmaHost = 'http://magma.test';
const projectName = 'example';
const rulesResponse: MagmaRulesResponse = exampleRulesResponse;
const existingNameCounter = 2;
const author = 'developer';
const isoTimeString = (new Date(1700499086877)).toISOString();
const existingName: MagmaListName = {
    identifier: `EXAMPLE-HS${existingNameCounter}`,
    author,
    name_created_at: isoTimeString,
};

interface MagmaBulkGenerateRequest {
    names: MagmaBulkGenerateName[]
}

const server = setupServer(
    rest.get(`${magmaHost}/gnomon/:projectName`, (req, res, ctx) => {
        return res(ctx.json(rulesResponse));
    }),

    rest.get(`${magmaHost}/gnomon/:projectName/list/subject`, (req, res, ctx) => {
        const existingNames: MagmaListName[] = [];

        if (req.url.searchParams.get('regex')?.indexOf(existingName.identifier) != -1) {
            existingNames.push(existingName);
        }

        return res(ctx.json(existingNames));
    }),

    rest.post(`${magmaHost}/gnomon/:projectName/generate`, (req, res, ctx) => {
        const payload = req.body as MagmaBulkGenerateRequest;

        return res(ctx.json({ created: payload.names, existing: [] }));
    }),
);

beforeAll(() => server.listen());
afterEach(() => {
    server.resetHandlers();
    jest.restoreAllMocks();
});
afterAll(() => server.close());

describe('NamesCreate', () => {
    it('fetches rules and can add names', async () => {

        const { renderResult } = renderWithProviders(
            <NamesCreate project_name={projectName} />
        );

        await waitFor(() => {
            expect(document.querySelector(
                'button[aria-label="Add Names"]:not(.Mui-disabled)'
            )).toBeTruthy();
        });
        expect(renderResult.asFragment()).toMatchSnapshot();

        const rule = 'biospecimen';

        // open menu; add a rule name
        fireEvent.click(screen.getByLabelText('Add Names'));
        fireEvent.click(screen.getByText(rule));

        // wait for menu to close; wait for rule group
        await waitForElementToBeRemoved(() => document.getElementById('add-names-rules-menu'));
        await screen.findByText(rule);

        expect(renderResult.asFragment()).toMatchSnapshot();
    });

    it('tracks name completion as elements are filled in', async () => {

        const { renderResult, store } = renderWithProviders(
            <NamesCreate project_name={projectName} />
        );

        await waitFor(() => {
            expect(document.querySelector(
                'button[aria-label="Add Names"]:not(.Mui-disabled)'
            )).toBeTruthy();
        });

        const rule = 'subject';

        // open menu; add a rule name
        fireEvent.click(screen.getByLabelText('Add Names'));
        fireEvent.click(screen.getByText(rule));

        // wait for menu to close; wait for rule group
        await waitForElementToBeRemoved(() => document.getElementById('add-names-rules-menu'));
        await screen.findByText(rule);

        expect(Object.keys(store.getState().names.completeCreateNames.byLocalId).length).toBe(0);

        fireEvent.click(screen.getByLabelText('Select Value for species'));
        fireEvent.click(screen.getByText('HS - Homo sapiens'));
        fireEvent.change(screen.getByPlaceholderText('n'), { target: { value: '1' } });

        expect(renderResult.asFragment()).toMatchSnapshot();
        expect(Object.keys(store.getState().names.completeCreateNames.byLocalId).length).toBe(1);

        await waitFor(() => {
            expect(document.querySelector(
                'div.hasError span[class*="infoTooltipIconContainer"]'
            )).toBeFalsy();

            const composeErrors = Object.entries(store.getState().names.createNameGroups.composeErrorsByLocalId);
            expect(composeErrors.length).toBe(1);
            expect(composeErrors[0][1].error).toBe(false);
        });
    });

    it('shows an error when the name is a duplicate', async () => {

        const { renderResult, store } = renderWithProviders(
            <NamesCreate project_name={projectName} />
        );

        await waitFor(() => {
            expect(document.querySelector(
                'button[aria-label="Add Names"]:not(.Mui-disabled)'
            )).toBeTruthy();
        });

        const rule = 'subject';

        // open menu; add a rule name
        fireEvent.click(screen.getByLabelText('Add Names'));
        fireEvent.click(screen.getByText(rule));

        // wait for menu to close; wait for rule group
        await waitForElementToBeRemoved(() => document.getElementById('add-names-rules-menu'));
        await screen.findByText(rule);

        expect(Object.keys(store.getState().names.completeCreateNames.byLocalId).length).toBe(0);

        fireEvent.click(screen.getByLabelText('Select Value for species'));
        fireEvent.click(screen.getByText('HS - Homo sapiens'));
        fireEvent.change(screen.getByPlaceholderText('n'), { target: { value: String(existingNameCounter) } });

        expect(renderResult.asFragment()).toMatchSnapshot();
        expect(Object.keys(store.getState().names.completeCreateNames.byLocalId).length).toBe(1);

        await waitFor(() => {
            expect(document.querySelector(
                'div.hasError span[class*="infoTooltipIconContainer"]'
            )).toBeTruthy();

            const composeErrors = Object.entries(store.getState().names.createNameGroups.composeErrorsByLocalId);
            expect(composeErrors.length).toBe(1);
            expect(composeErrors[0][1].error).toBe(true);
        });
    });

    it('selects names based on search criteria', async () => {

        /*
        1. create two names with different subject values
        2. do a search with one value
        3. assert only one name is selected
        */

        const { renderResult, store } = renderWithProviders(
            <NamesCreate project_name={projectName} />
        );

        await waitFor(() => {
            expect(document.querySelector(
                'button[aria-label="Add Names"]:not(.Mui-disabled)'
            )).toBeTruthy();
        });

        const rule = 'subject';

        // open menu; add a rule name
        fireEvent.click(screen.getByLabelText('Add Names'));
        fireEvent.click(screen.getByText(rule));

        // wait for menu to close; wait for rule group
        await waitForElementToBeRemoved(() => document.getElementById('add-names-rules-menu'));
        await screen.findByText(rule);

        fireEvent.click(screen.getByLabelText('Copy Name with Values'));

        await waitFor(() => {
            expect(document.querySelectorAll(
                'div[class*="composerRowContainer"]'
            ).length).toBe(2);
        });

        const selectButtons = screen.getAllByLabelText('Select Value for species');
        const counterInputs = screen.getAllByPlaceholderText('n');

        fireEvent.click(selectButtons[0]);
        fireEvent.click(screen.getByText('HS - Homo sapiens'));
        fireEvent.change(counterInputs[0], { target: { value: String(existingNameCounter + 1) } });

        await waitFor(() => {
            const composeErrors = Object.values(store.getState().names.createNameGroups.composeErrorsByLocalId);
            expect(composeErrors.length).toBe(2);
            expect(composeErrors.map(el => el.checkStatus)).toStrictEqual(['idle', 'idle']);
        });

        fireEvent.click(selectButtons[1]);
        fireEvent.click(screen.getByText('MM - Mus musculus'));
        fireEvent.change(counterInputs[1], { target: { value: String(existingNameCounter + 1) } });

        await waitFor(() => {
            const composeErrors = Object.values(store.getState().names.createNameGroups.composeErrorsByLocalId);
            expect(composeErrors.length).toBe(2);
            expect(composeErrors.map(el => el.checkStatus)).toStrictEqual(['idle', 'idle']);
        });

        // enter search criteria
        fireEvent.click(screen.getByLabelText('Find and Filter'));
        fireEvent.click(screen.getByLabelText('Select Value for Common Element'));
        // @ts-ignore
        fireEvent.click(document.querySelector(
            'li[role="menuitem"][aria-disabled="false"]')
        );

        await waitFor(() => {
            expect(document.querySelector(
                '#find-and-filter-dialogue div[class*="elementsEditorContainer"]'
            )).toBeTruthy();
        });

        // @ts-ignore
        fireEvent.click(document.querySelector(
            '#find-and-filter-dialogue [aria-label="Select Value for species"]'
        ));
        fireEvent.click(screen.getByText('HS - Homo sapiens'));

        // do search
        fireEvent.click(screen.getByText('Find + Select'));

        // wait for only one name to be selected
        await waitFor(() => {
            const checkboxStates = screen.getAllByLabelText('Select Name')
                // @ts-ignore
                .map(checkbox => checkbox.checked);

            expect(checkboxStates).toStrictEqual([true, false]);
        });
    });

    it('filters names based on search criteria', async () => {

        /*
        1. create two names with different subject values
        2. do a search with one value
        3. assert only one name is shown
        */

        const { renderResult, store } = renderWithProviders(
            <NamesCreate project_name={projectName} />
        );

        await waitFor(() => {
            expect(document.querySelector(
                'button[aria-label="Add Names"]:not(.Mui-disabled)'
            )).toBeTruthy();
        });

        const rule = 'subject';

        // open menu; add a rule name
        fireEvent.click(screen.getByLabelText('Add Names'));
        fireEvent.click(screen.getByText(rule));

        // wait for menu to close; wait for rule group
        await waitForElementToBeRemoved(() => document.getElementById('add-names-rules-menu'));
        await screen.findByText(rule);

        fireEvent.click(screen.getByLabelText('Copy Name with Values'));

        await waitFor(() => {
            expect(document.querySelectorAll(
                'div[class*="composerRowContainer"]'
            ).length).toBe(2);
        });

        const selectButtons = screen.getAllByLabelText('Select Value for species');
        const counterInputs = screen.getAllByPlaceholderText('n');

        fireEvent.click(selectButtons[0]);
        fireEvent.click(screen.getByText('HS - Homo sapiens'));
        fireEvent.change(counterInputs[0], { target: { value: String(existingNameCounter + 1) } });

        await waitFor(() => {
            const composeErrors = Object.values(store.getState().names.createNameGroups.composeErrorsByLocalId);
            expect(composeErrors.length).toBe(2);
            expect(composeErrors.map(el => el.checkStatus)).toStrictEqual(['idle', 'idle']);
        });

        fireEvent.click(selectButtons[1]);
        fireEvent.click(screen.getByText('MM - Mus musculus'));
        fireEvent.change(counterInputs[1], { target: { value: String(existingNameCounter + 1) } });

        await waitFor(() => {
            const composeErrors = Object.values(store.getState().names.createNameGroups.composeErrorsByLocalId);
            expect(composeErrors.length).toBe(2);
            expect(composeErrors.map(el => el.checkStatus)).toStrictEqual(['idle', 'idle']);
        });

        // enter search criteria
        fireEvent.click(screen.getByLabelText('Find and Filter'));
        fireEvent.click(screen.getByLabelText('Select Value for Common Element'));
        // @ts-ignore
        fireEvent.click(document.querySelector(
            'li[role="menuitem"][aria-disabled="false"]')
        );

        await waitFor(() => {
            expect(document.querySelector(
                '#find-and-filter-dialogue div[class*="elementsEditorContainer"]'
            )).toBeTruthy();
        });

        // @ts-ignore
        fireEvent.click(document.querySelector(
            '#find-and-filter-dialogue [aria-label="Select Value for species"]'
        ));
        fireEvent.click(screen.getByText('HS - Homo sapiens'));

        // do search
        fireEvent.click(screen.getByText('Filter'));

        // wait for only one name to be in the DOM
        await waitFor(() => {
            expect(document.querySelectorAll(
                'div[class*="composerRowContainer"]'
            ).length).toBe(1);
        });
    });

    it('can copy names n times', async () => {

        const { renderResult, store } = renderWithProviders(
            <NamesCreate project_name={projectName} />
        );

        await waitFor(() => {
            expect(document.querySelector(
                'button[aria-label="Add Names"]:not(.Mui-disabled)'
            )).toBeTruthy();
        });

        const rule = 'subject';

        // open menu; add a rule name
        fireEvent.click(screen.getByLabelText('Add Names'));
        fireEvent.click(screen.getByText(rule));

        // wait for menu to close; wait for rule group
        await waitForElementToBeRemoved(() => document.getElementById('add-names-rules-menu'));
        await screen.findByText(rule);

        fireEvent.click(screen.getByLabelText('Select Value for species'));
        fireEvent.click(screen.getByText('HS - Homo sapiens'));

        expect(renderResult.asFragment()).toMatchSnapshot();
        expect(Object.keys(store.getState().names.completeCreateNames.byLocalId).length).toBe(0);

        // select the name
        fireEvent.click(screen.getByLabelText('Select Name'));
        await waitFor(() => {
            // @ts-ignore
            expect(screen.getByLabelText('Select Name').checked).toBeTruthy();
            // @ts-ignore
            expect(screen.getByLabelText('Copy and Replace').disabled).toBeFalsy();
        });

        fireEvent.click(screen.getByLabelText('Copy and Replace'));

        const copyQuantity = 4;
        fireEvent.change(screen.getByLabelText('copy-radio-quantity'), { target: { value: String(copyQuantity) } });
        fireEvent.click(screen.getByLabelText('perform-copy-or-replace'));

        await waitFor(() => {
            expect(document.querySelectorAll(
                'div[class*="composerRowContainer"]'
            ).length).toBe(1 + copyQuantity);

            expect(
                screen.queryAllByText('HS').length
            ).toBe(1 + copyQuantity);
        });
    });

    it('can iterate on a rule and replace values', async () => {

        const { renderResult, store } = renderWithProviders(
            <NamesCreate project_name={projectName} />
        );

        await waitFor(() => {
            expect(document.querySelector(
                'button[aria-label="Add Names"]:not(.Mui-disabled)'
            )).toBeTruthy();
        });

        const rule = 'subject';
        const firstNameCounterValue = existingNameCounter + 1;
        const iterateDelta = 5;

        // open menu; add a rule name
        fireEvent.click(screen.getByLabelText('Add Names'));
        fireEvent.click(screen.getByText(rule));

        // wait for menu to close; wait for rule group
        await waitForElementToBeRemoved(() => document.getElementById('add-names-rules-menu'));
        await screen.findByText(rule);

        fireEvent.click(screen.getByLabelText('Select Value for species'));
        fireEvent.click(screen.getByText('HS - Homo sapiens'));
        fireEvent.change(screen.getByPlaceholderText('n'), { target: { value: String(firstNameCounterValue) } });

        expect(renderResult.asFragment()).toMatchSnapshot();
        expect(Object.keys(store.getState().names.completeCreateNames.byLocalId).length).toBe(1);

        // select the name
        fireEvent.click(screen.getByLabelText('Select Name'));
        await waitFor(() => {
            // @ts-ignore
            expect(screen.getByLabelText('Select Name').checked).toBeTruthy();
            // @ts-ignore
            expect(screen.getByLabelText('Copy and Replace').disabled).toBeFalsy();
        });

        fireEvent.click(screen.getByLabelText('Copy and Replace'));
        fireEvent.click(screen.getByText('Iterate on Element'));
        fireEvent.click(screen.getByLabelText('Select Value for Element'));
        // click the 'subject' menu item
        // @ts-ignore
        fireEvent.click(document.querySelector(
            'li[role="menuitem"][aria-disabled="false"]')
        );
        fireEvent.change(
            screen.getByLabelText('iterate-radio-boundary-minimum'),
            { target: { value: String(firstNameCounterValue + 1) } },
        );
        fireEvent.change(
            screen.getByLabelText('iterate-radio-boundary-maximum'),
            { target: { value: String(firstNameCounterValue + 1 + iterateDelta) } },
        );

        fireEvent.click(screen.getByLabelText('perform-copy-or-replace'));

        await waitFor(() => {
            const expectedNamesCount = iterateDelta + 2;

            const composeErrors = Object.values(store.getState().names.createNameGroups.composeErrorsByLocalId);
            expect(composeErrors.length).toBe(expectedNamesCount);
            expect(_.every(composeErrors.map(el => el.error))).toBe(false);

            const nameComposers = document.querySelectorAll(
                'div[class*="composerRowContainer"]'
            );
            expect(nameComposers.length).toBe(expectedNamesCount);

            nameComposers.forEach((nameComposer, idx) => {
                // @ts-ignore
                expect(nameComposer.querySelector(
                    '[aria-label="Select Value for species"] .MuiButton-label'
                ).textContent).toBe('HS');

                // @ts-ignore
                expect(nameComposer.querySelector(
                    '[aria-label="subject-counter-value"]'
                    // @ts-ignore
                ).value).toBe(String(idx + firstNameCounterValue));
            });
        });

        // replace values
        fireEvent.click(screen.getByLabelText('Select the Name Group'));
        fireEvent.click(screen.getByLabelText('Copy and Replace'));
        fireEvent.click(screen.getByText('Replace Values'));
        // @ts-ignore
        fireEvent.click(document.querySelector(
            '#copy-and-replace-dialogue [aria-label="Select Value for species"]'
        ));
        fireEvent.click(screen.getByText('MM - Mus musculus'));

        fireEvent.click(screen.getByLabelText('perform-copy-or-replace'));

        await waitFor(() => {
            const expectedNamesCount = iterateDelta + 2;

            const composeErrors = Object.values(store.getState().names.createNameGroups.composeErrorsByLocalId);
            expect(composeErrors.length).toBe(expectedNamesCount);
            expect(_.every(composeErrors.map(el => el.error))).toBe(false);

            const nameComposers = document.querySelectorAll(
                'div[class*="composerRowContainer"]'
            );
            expect(nameComposers.length).toBe(expectedNamesCount);

            nameComposers.forEach((nameComposer, idx) => {
                // @ts-ignore
                expect(nameComposer.querySelector(
                    '[aria-label="Select Value for species"] .MuiButton-label'
                ).textContent).toBe('MM');

                // @ts-ignore
                expect(nameComposer.querySelector(
                    '[aria-label="subject-counter-value"]'
                    // @ts-ignore
                ).value).toBe(String(idx + firstNameCounterValue));
            });
        });
    });

    it('can bulk delete', async () => {

        const { renderResult, store } = renderWithProviders(
            <NamesCreate project_name={projectName} />
        );

        await waitFor(() => {
            expect(document.querySelector(
                'button[aria-label="Add Names"]:not(.Mui-disabled)'
            )).toBeTruthy();
        });

        const rule = 'subject';

        // open menu; add a rule name
        fireEvent.click(screen.getByLabelText('Add Names'));
        fireEvent.click(screen.getByText(rule));

        // wait for menu to close; wait for rule group
        await waitForElementToBeRemoved(() => document.getElementById('add-names-rules-menu'));
        await screen.findByText(rule);

        fireEvent.click(screen.getByLabelText('Copy Name with Values'));

        await waitFor(() => {
            expect(document.querySelectorAll(
                'div[class*="composerRowContainer"]'
            ).length).toBe(2);
        });

        // select names
        fireEvent.click(screen.getByLabelText('Select the Name Group'));
        await waitFor(() => {
            expect(
                // @ts-ignore
                _.every(screen.getAllByLabelText('Select Name').map(el => el.checked))
            ).toBe(true);
            // @ts-ignore
            expect(screen.getByLabelText('Delete Selection').disabled).toBeFalsy();
        });


        // do deletion
        fireEvent.click(screen.getByLabelText('Delete Selection'));
        fireEvent.click(screen.getByLabelText('confirm-confirmation'));

        await waitFor(() => {
            expect(document.querySelectorAll(
                'div[class*="composerRowContainer"]'
            ).length).toBe(0);
        });
    });

    it('creates complete names and exports to file', async () => {

        const { renderResult, store } = renderWithProviders(
            <NamesCreate project_name={projectName} />
        );

        await waitFor(() => {
            expect(document.querySelector(
                'button[aria-label="Add Names"]:not(.Mui-disabled)'
            )).toBeTruthy();
        });

        const rule = 'subject';

        // open menu; add a rule name
        fireEvent.click(screen.getByLabelText('Add Names'));
        fireEvent.click(screen.getByText(rule));

        // wait for menu to close; wait for rule group
        await waitForElementToBeRemoved(() => document.getElementById('add-names-rules-menu'));
        await screen.findByText(rule);

        fireEvent.click(screen.getByLabelText('Copy Name with Values'));

        await waitFor(() => {
            expect(document.querySelectorAll(
                'div[class*="composerRowContainer"]'
            ).length).toBe(2);
        });

        const selectButtons = screen.getAllByLabelText('Select Value for species');
        const counterInputs = screen.getAllByPlaceholderText('n');

        fireEvent.click(selectButtons[0]);
        fireEvent.click(screen.getByText('HS - Homo sapiens'));
        fireEvent.change(counterInputs[0], { target: { value: String(existingNameCounter + 1) } });

        await waitFor(() => {
            const composeErrors = Object.values(store.getState().names.createNameGroups.composeErrorsByLocalId);
            expect(composeErrors.length).toBe(2);
            expect(_.every(composeErrors.map(el => el.error))).toBe(false);
        });

        fireEvent.click(selectButtons[1]);
        fireEvent.click(screen.getByText('MM - Mus musculus'));
        fireEvent.change(counterInputs[1], { target: { value: String(existingNameCounter + 1) } });

        await waitFor(() => {
            const composeErrors = Object.values(store.getState().names.createNameGroups.composeErrorsByLocalId);
            expect(composeErrors.length).toBe(2);
            expect(_.every(composeErrors.map(el => el.error))).toBe(false);
        });


        interface Row {
            name: string
            rule: string
        }

        const completedNames: Row[] = selectCompleteCreateNamesCreationPayloads(store.getState())
            .map(payload => ({ name: payload.renderedName, rule: payload.ruleName }));

        // do creation
        fireEvent.click(screen.getByLabelText('Create'));

        await waitFor(() => {
            for (const completedName of completedNames) {
                expect(screen.getByText(completedName.name)).toBeTruthy();
            }
        });

        fireEvent.click(screen.getByLabelText('perform-creation'));

        await waitFor(() => {
            const creationRequestState = store.getState().names.creationRequest;

            expect(creationRequestState.status).toBe('success');
            expect(screen.getByText('Export xlsx')).toBeTruthy();
        });

        // do export
        const exportSpy = jest.spyOn(exportModule, 'exportDataToBlob');
        const createObjectURLMock = jest.fn();
        window.URL.createObjectURL = createObjectURLMock;

        const expectedBlob = await exportModule.exportDataToBlob(completedNames, 'xlsx');

        fireEvent.click(screen.getByText('Export xlsx'));

        await waitFor(() => expect(exportSpy).toHaveBeenCalledWith(completedNames, 'xlsx'));
        await waitFor(() => expect(createObjectURLMock).toHaveBeenCalledWith(expectedBlob));

        await waitFor(() => {
            expect(document.querySelector(
                'span.success[class*="requestStatus"]'
            )).toBeTruthy();
        });
    });
});