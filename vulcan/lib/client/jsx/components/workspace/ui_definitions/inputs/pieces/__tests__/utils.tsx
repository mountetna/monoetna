import React from 'react';
import { render, RenderResult } from '@testing-library/react';
import { ThemeProvider } from '@material-ui/core/styles';

import {createEtnaTheme} from 'etna-js/style/theme';
const theme = createEtnaTheme('#ffaa44', '#948f8e');

interface Render {
    renderResult: RenderResult
}

export function renderWithProviders(
    component: React.ReactElement
): Render {

    const renderResult = render(
        <ThemeProvider theme={theme}>
            {component}
        </ThemeProvider>
    );

    return { renderResult };
}