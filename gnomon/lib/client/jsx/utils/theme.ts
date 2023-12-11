import { createTheme, Theme } from '@material-ui/core/styles';
// @ts-ignore
import Color from 'color';



export function createCustomColorButtonTheme(color: string, baseTheme?: Theme): Theme {
    const primary = Color(color);

    return createTheme(baseTheme, {
        palette: {
            primary: {
                main: primary.string()
            }
        },
        overrides: {
            MuiButton: {
                containedPrimary: {
                    backgroundColor: primary.alpha(0.3).string(),
                    color: primary.darken(0.5).string(),
                    '&:hover': {
                        backgroundColor: primary.alpha(0.5).string()
                    },
                    '&:disabled': {
                        backgroundColor: primary.alpha(0.15).string()
                    },
                },
            },
        }
    });
}