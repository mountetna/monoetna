import { ThemeProvider, createMuiTheme } from '@material-ui/core/styles';
import Color from 'color';

export const createEtnaTheme = (p,s) => {
  const primary = Color(p);
  const secondary = Color(s);
  const theme = {
    typography: {
      fontFamily: 'Open Sans,sans-serif'
    },
    shape: {
      borderRadius: "2px"
    },
    palette: {
      primary: {
        main: primary.string()
      },
      secondary: {
        main: secondary.string()
      }
    },
    overrides: {
      MuiButton: {
        root: {
          textTransform: "none"
        },
        containedPrimary: {
          backgroundColor: primary.alpha(0.3).string(),
          color: primary.darken(0.5).string(),
          '&:hover': {
            backgroundColor: primary.alpha(0.5).string()
          }
        },
        containedSecondary: {
          backgroundColor: secondary.alpha(0.3).string(),
          color: secondary.darken(0.5).string(),
          '&:hover': {
            backgroundColor: secondary.alpha(0.5).string()
          }
        },
      },
      MuiChip: {
        root: {
          borderRadius: "4px",
          boxShadow: "0 0 3px #ccc",
          margin: "0px 2px"
        }
      },
      MuiTableCell: {
        root: {
          fontSize: "1rem"
        }
      },
      MuiIconButton: {
        root: {
          borderRadius: "2px"
        }
      },
      MuiScopedCssBaseline: {
        root: {
          backgroundColor: "none"
        }
      },
      MuiTableRow: {
        head: {
          fontFamily: 'Cousine, monospace',
          color: '#333',
          background: '#eee'
        }
      }
    },
    props: {
      MuiButton: {
        size: "small",
        variant: "contained",
        color: "primary",
        disableElevation: true,
        disableRipple: true
      },
      MuiCheckbox: {
        disableRipple: true
      },
      MuiChip: {
        color: "secondary",
        size: "small"
      },
      MuiTableContainer: {
        disableElevation: true,
        variant: "outlined"
      },
      MuiPaper: {
        square: true
      },
      MuiTable: {
        size: "small"
      },
      MuiAutocomplete: {
        size: "small"
      }
    }
  };
  return createMuiTheme(theme);
}
