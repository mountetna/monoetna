import React from 'react';
import { rest } from 'msw';
import { setupServer } from 'msw/node';
import { fireEvent, waitFor, waitForElementToBeRemoved, screen } from '@testing-library/react';

import NamesCreate from '../names-create';
import { renderWithProviders } from '../../../utils/test';
import * as exampleRulesResponse from './static/example-rules-response.json';
import { MagmaRulesResponse } from '../../../utils/rules';
import { MagmaListName } from '../../../utils/names';


global.ResizeObserver = jest.fn().mockImplementation(() => ({
    observe: jest.fn(),
    unobserve: jest.fn(),
    disconnect: jest.fn(),
}));

const magmaHost = 'http://magma.test';
const projectName = 'example';
const rulesResponse: MagmaRulesResponse = exampleRulesResponse;
const existingNameCounter = 2;
const existingName: MagmaListName = {
    identifier: `EXAMPLE-HS${existingNameCounter}`,
    author: 'developer',
    name_created_at: (new Date(1700499086877)).toISOString(),
};

const server = setupServer(
    rest.get(`${magmaHost}/gnomon/:projectName`, (req, res, ctx) => {
        return res(ctx.json(rulesResponse));
    }),

    rest.get(`${magmaHost}/gnomon/:projectName/list/subject`, (req, res, ctx) => {
        const url = new URL(req.url);
        const existingNames: MagmaListName[] = [];

        if (url.searchParams.get('regex')?.indexOf(existingName.identifier) != -1) {
            existingNames.push(existingName);
        }

        return res(ctx.json(existingNames));
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
        // await waitFor(() => {
        //     expect(screen.queryByText(rule)).toBeTruthy();
        // });
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
        // await waitFor(() => {
        //     expect(screen.queryByText(rule)).toBeTruthy();
        // });
        fireEvent.click(screen.getByText(rule));

        // wait for menu to close; wait for rule group
        await waitForElementToBeRemoved(() => document.getElementById('add-names-rules-menu'));
        await screen.findByText(rule);

        expect(Object.keys(store.getState().names.completeCreateNames.byLocalId).length == 0);

        fireEvent.click(screen.getByLabelText('Select Value for species'));
        fireEvent.click(screen.getByText('HS - Homo sapiens'));
        fireEvent.change(screen.getByPlaceholderText('n'), { target: { value: '1' } });

        expect(renderResult.asFragment()).toMatchSnapshot();
        expect(Object.keys(store.getState().names.completeCreateNames.byLocalId).length == 1);

        await waitFor(() => {
            expect(document.querySelector(
                'div.hasError span[class*="infoTooltipIconContainer"]'
            )).toBeFalsy();
            
            const composeErrors = Object.entries(store.getState().names.createNameGroups.composeErrorsByLocalId);
            expect(composeErrors.length).toBe(1);
            expect(composeErrors[0][1]).toBe(false);
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
            expect(composeErrors[0][1]).toBe(true);
        });
    });
});