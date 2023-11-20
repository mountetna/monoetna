import React from 'react'
import { Provider } from 'react-redux';
import { render, RenderResult } from '@testing-library/react';
import { ThemeProvider } from '@material-ui/core/styles';

import createStore, { State } from '../store';
import { theme } from '../gnomon-ui';



export function renderWithProviders(
    component: React.ReactElement,
    preloadedState?: Partial<State>
): RenderResult {

    return render(
        <Provider store={createStore(preloadedState)}>
            <ThemeProvider theme={theme}>
                {component}
            </ThemeProvider>
        </Provider>
    )
}