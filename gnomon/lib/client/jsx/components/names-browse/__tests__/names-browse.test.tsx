import React from 'react';
import { rest } from 'msw';
import { setupServer } from 'msw/node';
import { fireEvent, waitFor, waitForElementToBeRemoved, screen } from '@testing-library/react';

import NamesBrowse, { Name } from '../names-browse';
import { renderWithProviders } from '../../../utils/test';
import { MagmaRules } from '../../../utils/rules';
import { MagmaListName } from '../../../utils/names';
import * as exportModule from '../../../utils/export';



const magmaHost = 'http://magma.test';
const projectName = 'test_project';
const rules = ['rule1', 'rule2'];
const namesByRule: Record<string, MagmaListName[]> = {};

rules.forEach((rule, idx) => {
    const names: MagmaListName[] = [`name${idx * 2 + 1}`, `name${idx * 2 + 2}`].map(name => ({
        identifier: name,
        author: 'test_author',
        name_created_at: (new Date(1700499086877)).toISOString(),
    }));
    namesByRule[rule] = names;
});

const server = setupServer(
    rest.get(`${magmaHost}/gnomon/${projectName}`, (req, res, ctx) => {
        const _rules: Record<string, string> = {};
        Object.keys(namesByRule).forEach(rule => _rules[rule] = 'TOKEN');

        const rules: MagmaRules = {
            rules: _rules,
            tokens: { TOKEN: { name: 'TOKEN', label: 'token', values: { VALUE: 'value' } } },
            synonyms: [],
        };
        return res(ctx.json({ config: rules }));
    }),

    ...Object.entries(namesByRule).map(([rule, names]) => (
        rest.get(`${magmaHost}/gnomon/${projectName}/list/${rule}`, (req, res, ctx) => {
            return res(ctx.json(names));
        })
    )),
);

beforeAll(() => server.listen());
afterEach(() => {
    server.resetHandlers();
    jest.restoreAllMocks();
});
afterAll(() => server.close());


describe('NamesBrowse', () => {
    it('fetches and renders names for all rules', async () => {

        const { renderResult } = renderWithProviders(
            <NamesBrowse project_name={projectName} />
        );

        await waitFor(() => screen.getByText(/name4/));
        expect(renderResult.asFragment()).toMatchSnapshot();
    });

    it('exports all names to csv', async () => {

        const { renderResult } = renderWithProviders(
            <NamesBrowse project_name={projectName} />
        );

        await waitFor(() => screen.findByText(/name4/));

        const exportSpy = jest.spyOn(exportModule, 'exportDataToBlob');
        const createObjectURLMock = jest.fn();
        window.URL.createObjectURL = createObjectURLMock;

        const expectedData: Name[] = [];
        Object.entries(namesByRule).forEach(([rule, names]) => {
            names.forEach(name => {
                expectedData.push({
                    name: name.identifier,
                    author: name.author,
                    ruleName: rule,
                    createdAt: name.name_created_at,
                });
            });
        });
        const expectedBlob = await exportModule.exportDataToBlob(expectedData, 'csv');

        fireEvent.click(screen.getByLabelText('Export'));
        fireEvent.click(screen.getByText('csv'));

        await waitForElementToBeRemoved(() => document.getElementById('export-file-formats'));
        await waitFor(() => expect(exportSpy).toHaveBeenCalledWith(expectedData, 'csv'));
        await waitFor(() => expect(createObjectURLMock).toHaveBeenCalledWith(expectedBlob));
    });
});
