import * as React from 'react';
import { ThemeProvider, createMuiTheme } from '@material-ui/core/styles';

import JanusNav from './janus-nav';
import JanusMain from './janus-main';
import JanusAdmin from './janus-admin';
import JanusSettings from './janus-settings';
import ProjectView from './project-view';
import FlagsView from './flags/flags-view';

import { findRoute, setRoutes } from 'etna-js/dispatchers/router';

const theme = createMuiTheme({
  typography: {
    fontFamily: 'Open Sans,sans-serif'
  },
  shape: {
    borderRadius: "2px"
  },
  palette: {
    primary: {
      main: "#3684fd"
    },
    secondary: {
      main: "#77c"
    }
  },
  overrides: {
    MuiButton: {
      root: {
        textTransform: "none"
      }
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
    MuiGrid: {
      disableElevation: true
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
    }
  }
});

const ROUTES = [
  {
    name: 'main',
    template: '',
    component: JanusMain
  },
  {
    name: 'admin',
    template: 'admin',
    component: JanusAdmin
  },
  {
    name: 'settings',
    template: 'settings',
    component: JanusSettings
  },
  {
    name: 'flags',
    template: 'flags',
    component: FlagsView
  },
  {
    name: 'projects',
    template: ':project_name',
    component: ProjectView
  }
];

setRoutes(ROUTES);

const Invalid = () => <div>Path invalid</div>;

const JanusUI = () => {
  let { route, params }  = findRoute({ path: window.location.pathname }, ROUTES);
  let Component = route ? route.component : Invalid;

  return (
    <ThemeProvider theme={theme}>
      <div id='janus-group'>
        <JanusNav/>
        <Component {...params}/>
      </div>
    </ThemeProvider>
  );
}

export default JanusUI;
