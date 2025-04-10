import React from 'react';
import { fireEvent, waitFor, screen } from '@testing-library/react';
import { renderWithProviders } from './utils';
import { ComponentUse } from '../../components/visualizations';
import { nestedOptionSet, keyedValues } from '../utils';
import NestedDropdownMultiChoiceAdvancedPiece from '../nested_dropdown_multi_choice_advanced_piece';

const initialValues: keyedValues = {
    marks_use: ['cd3'],
    some_other: 3
};
let values = initialValues;

const optionSets: nestedOptionSet = {
    Row_Features: {
        set1: {
            a: null,
            b: null,
            cd400: null,
            'cd4.1': null
        },
        set2: {
            cd3: null,
            cd4: null
        }
    },
    Col_Features: {
        d: null,
        e: null,
        f: null
    }
};
const updateValue = (newValue: any, key: string, prevValues = {...values}) => {
    prevValues[key] = newValue;
    values = prevValues;
};

describe('NestedDropdownMultiChoicePieces', () => {
    beforeEach(() => {
        values = initialValues;
    });
    describe('without additions', () => {
        beforeEach(async () => {
            const { renderResult } = renderWithProviders(
                <ComponentUse
                    key='marks_use'
                    k='marks_use'
                    value={values.marks_use}
                    extra_inputs={[
                        'Primary Data',
                        optionSets,
                        false,
                        false
                    ]}
                    updateValue={updateValue}
                    comps={{marks_use: NestedDropdownMultiChoiceAdvancedPiece}}
                />
            );
            await waitFor(() => {
                expect(document.querySelector(
                    'input[id="MultiMultiChoice-optionPaths"]:not(.Mui-disabled)'
                )).toBeTruthy();
            });
        });

        it('shows chosen values', async () => {
            expect(values.marks_use).toEqual(initialValues.marks_use);
            // Dropdown for matching optionSet exists
            expect(document.querySelector('input[id="MultiMultiChoice-Row_Features---set2-leaves"]:not(.Mui-disabled)')).toBeTruthy();
            // values' chips exist
            expect(screen.getByText('cd3')).toBeTruthy();
        });

        it('does not give bulk add or reorder', () => {
            expect(screen.queryByTitle('Bulk Add')).toBeFalsy();
            expect(screen.queryByText('Reorder selections?')).toBeFalsy();
        });
    });

    describe('with Bulk Add', () => {
        beforeEach(async () => {
            const { renderResult } = renderWithProviders(
                <ComponentUse
                    key='marks_use'
                    k='marks_use'
                    value={values.marks_use}
                    extra_inputs={[
                        'Primary Data',
                        optionSets,
                        false,
                        true
                    ]}
                    updateValue={updateValue}
                    comps={{marks_use: NestedDropdownMultiChoiceAdvancedPiece}}
                />
            );
            await waitFor(() => {
                expect(document.querySelector(
                    'input[id="MultiMultiChoice-optionPaths"]:not(.Mui-disabled)'
                )).toBeTruthy();
            });
        });

        it('shows chosen values', async () => {
            expect(values.marks_use).toEqual(initialValues.marks_use);
            // Dropdown for matching optionSet exists
            expect(document.querySelector('input[id="MultiMultiChoice-Row_Features---set2-leaves"]:not(.Mui-disabled)')).toBeTruthy();
            // values' chips exist
            expect(screen.getByText('cd3')).toBeTruthy();
        });

        it('shows bulk add', () => {
            expect(screen.queryByTitle('Bulk Add')).toBeTruthy();
            expect(screen.queryByText('Reorder selections?')).toBeFalsy();
        });

        it('merges matches from Bulk Add', async () => {
            // Open Bulk Add Modal, searches on given text input, accept default match selections
            fireEvent.click(screen.getByTitle('Bulk Add'));
            await waitFor(() => document.getElementById('TextField-input'));
            fireEvent.change(screen.getByTestId('bulk-add-user-text-input'), {target: {value: 'a,'}});
            fireEvent.click(screen.getByText('Search'));
            await waitFor(() => {
                expect(document.querySelector(
                    'button[aria-label="Add Selected Matches"]:not(.Mui-disabled)'
                )).toBeTruthy();
            });
            fireEvent.click(screen.getByText('Add Selected Matches'));

            // merges values
            expect(values.marks_use).toEqual(['cd3', 'a']);
        });

        it('allows user adjust number of matches given', async () => {
            // Open Bulk Add Modal, searches on given text input
            fireEvent.click(screen.getByTitle('Bulk Add'));
            await waitFor(() => document.getElementById('TextField-input'));
            fireEvent.change(screen.getByTestId('bulk-add-user-text-input'), {target: {value: 'cd4,'}});
            fireEvent.click(screen.getByText('Search'));
            await waitFor(() => {
                expect(document.querySelector(
                    'button[aria-label="Add Selected Matches"]:not(.Mui-disabled)'
                )).toBeTruthy();
            });
            expect(screen.getByText('cd4.1', {selector: 'span'})).toBeTruthy();
            expect(screen.getByText('cd400', {selector: 'span'})).toBeTruthy();

            // Decrease num matches to 1, check that 2nd match goes away
            fireEvent.click(screen.getByRole('slider', {name: 'Matches to Show'}));
            fireEvent.keyDown(screen.getByRole('slider', {name: 'Matches to Show'}), {key: 'ArrowLeft'});
            fireEvent.keyDown(screen.getByRole('slider', {name: 'Matches to Show'}), {key: 'ArrowLeft'});
            fireEvent.keyDown(screen.getByRole('slider', {name: 'Matches to Show'}), {key: 'ArrowLeft'});
            fireEvent.keyDown(screen.getByRole('slider', {name: 'Matches to Show'}), {key: 'ArrowLeft'});
            fireEvent.click(screen.getByText('Search'));
            await waitFor(() => {
                expect(document.querySelector(
                    'button[aria-label="Add Selected Matches"]:not(.Mui-disabled)'
                )).toBeTruthy();
            });
            expect(screen.getByText('cd4.1', {selector: 'span'})).toBeTruthy();
            expect(screen.queryByText('cd400', {selector: 'span'})).toBeFalsy();
        });

        it('allows user to select additional matches', async () => {
            // Open Bulk Add Modal, searches on given text input, accept default match selections
            fireEvent.click(screen.getByTitle('Bulk Add'));
            await waitFor(() => document.getElementById('TextField-input'));
            fireEvent.change(screen.getByTestId('bulk-add-user-text-input'), {target: {value: 'cd4,'}});
            fireEvent.click(screen.getByText('Search'));
            await waitFor(() => {
                expect(document.querySelector(
                    'button[aria-label="Add Selected Matches"]:not(.Mui-disabled)'
                )).toBeTruthy();
            });
            expect(screen.getByText('cd4', {selector: 'span'})).toBeTruthy();
            expect(screen.getByText('cd4.1', {selector: 'span'})).toBeTruthy();

            fireEvent.click(screen.getByText('cd4.1'));
            fireEvent.click(screen.getByText('Add Selected Matches'));

            // merges values
            expect(values.marks_use).toEqual(['cd3', 'cd4', 'cd4.1']);
        });
    });

    describe('with Reorder', () => {
        beforeEach(async () => {
            const { renderResult } = renderWithProviders(
                <ComponentUse
                    key='marks_use'
                    k='marks_use'
                    value={values.marks_use}
                    extra_inputs={[
                        'Primary Data',
                        optionSets,
                        true,
                        false
                    ]}
                    updateValue={updateValue}
                    comps={{marks_use: NestedDropdownMultiChoiceAdvancedPiece}}
                />
            );
            await waitFor(() => {
                expect(document.querySelector(
                    'input[id="MultiMultiChoice-optionPaths"]:not(.Mui-disabled)'
                )).toBeTruthy();
            });
        });

        it('shows chosen values', async () => {
            expect(values.marks_use).toEqual(initialValues.marks_use);
            // Dropdown for matching optionSet exists
            expect(document.querySelector('input[id="MultiMultiChoice-Row_Features---set2-leaves"]:not(.Mui-disabled)')).toBeTruthy();
            // values' chips exist
            expect(screen.getByText('cd3')).toBeTruthy();
        });

        it('shows reorder piece', () => {
            expect(screen.queryByTitle('Bulk Add')).toBeFalsy();
            expect(screen.queryByText('Reorder selections?')).toBeTruthy();
        });
    });
});

