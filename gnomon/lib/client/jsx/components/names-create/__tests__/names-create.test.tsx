import React from 'react';
import { rest } from 'msw'
import { setupServer } from 'msw/node'
import { fireEvent, waitFor, waitForElementToBeRemoved, screen } from '@testing-library/react';

import NamesCreate from '../names-create';
import { renderWithProviders } from '../../../utils/test';
import * as exampleRulesResponse from './static/example-rules-response.json';
import { MagmaRulesResponse } from '../../../utils/rules';


global.ResizeObserver = jest.fn().mockImplementation(() => ({
    observe: jest.fn(),
    unobserve: jest.fn(),
    disconnect: jest.fn(),
}))

const magmaHost = 'http://magma.test'
const projectName = 'example'
const rulesResponse: MagmaRulesResponse = exampleRulesResponse

const server = setupServer(
    rest.get(`${magmaHost}/gnomon/${projectName}`, (req, res, ctx) => {
        return res(ctx.json(rulesResponse));
    }),
)

beforeAll(() => server.listen());
afterEach(() => {
    server.resetHandlers();
    jest.restoreAllMocks();
});
afterAll(() => server.close());

describe('NamesCreate', () => {
    it("fetches rules and let's you add names", async () => {

        const component = renderWithProviders(
            <NamesCreate project_name={projectName} />
        )

        await waitFor(() => {
            expect(document.querySelector(
                'button[aria-label="Add Names"]:not(.Mui-disabled)'
            )).toBeTruthy()
        })
        expect(component.asFragment()).toMatchSnapshot()

        const rule = 'biospecimen'

        // open menu; add a rule name
        fireEvent.click(screen.getByLabelText('Add Names'))
        await waitFor(() => {
            expect(screen.queryByText(rule)).toBeTruthy();
        })
        fireEvent.click(screen.getByText(rule))

        // // // wait for menu to close; wait for rule group
        await waitForElementToBeRemoved(() => document.getElementById('add-names-rules-menu'))
        await screen.findByText(rule)

        expect(component.asFragment()).toMatchSnapshot()
    })

    
})