import React from 'react';
import { Provider } from 'react-redux';
import { Store } from 'redux';
import { render, RenderResult } from '@testing-library/react';
import { ThemeProvider } from '@material-ui/core/styles';

import createStore, { State } from '../store';
import { theme } from '../gnomon-ui';



interface Render {
    renderResult: RenderResult
    store: Store<State>
}



export function renderWithProviders(
    component: React.ReactElement,
    preloadedState?: Partial<State>
): Render {

    const store = createStore(preloadedState);

    const renderResult = render(
        <Provider store={store}>
            <ThemeProvider theme={theme}>
                {component}
            </ThemeProvider>
        </Provider>
    );

    return { renderResult, store };
}