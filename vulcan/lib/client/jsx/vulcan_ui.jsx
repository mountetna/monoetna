import React from 'react';
import {connect} from 'react-redux';
import {findRoute, setRoutes} from './router';

import {VulcanProvider} from './contexts/vulcan_context';
import { ThemeProvider } from '@material-ui/core/styles';

import { createEtnaTheme } from 'etna-js/style/theme';

const theme = createEtnaTheme("#de5833","#948f8e");

// Components.
import Browser from './components/browser.tsx';
import FigureList from './components/figure_list';
import VulcanNav from './components/vulcan_nav';
import Messages from 'etna-js/components/messages';
import {selectUser} from 'etna-js/selectors/user-selector';

import {showMessages} from 'etna-js/actions/message_actions';
import {updateLocation} from 'etna-js/actions/location_actions';
import RootView from 'etna-js/components/RootView';

import {ModalDialogContainer} from 'etna-js/components/ModalDialogContainer';
import {Notifications} from 'etna-js/components/Notifications';
import ReactModal from "react-modal";

const ROUTES = [
  {
    template: '',
    component: RootView,
    mode: ''
  },
  {
    template: ':project_name/',
    component: FigureList,
    mode: ''
  },
  {
    name: 'workflows',
    template: ':project_name/workflow',
    component: Browser,
    mode: 'workflow'
  },
  {
    name: 'workflow',
    template: ':projectName/workflow/:workflowName',
    component: Browser,
    mode: 'workflow'
  }
];

setRoutes(ROUTES);

const Empty = () => <div />;

class VulcanUI extends React.Component {
  constructor(props) {
    super(props);

    window.onpopstate = this.updateLocation.bind(this);
  }

  updateLocation() {
    let {updateLocation} = this.props;
    updateLocation(location);
  }

  componentDidMount() {
    ReactModal.setAppElement('#root');
  }

  render() {
    let {location, showMessages, environment, user} = this.props;
    let {route, params} = findRoute(location, ROUTES);
    let Component;
    let mode;

    if (this.props.Component) {
      mode = 'home';
      Component = this.props.Component;
    } else if (!route) {
      showMessages(['### You have lost your way: Path Invalid.']);
      mode = 'home';
      Component = Empty;
    } else {
      mode = route.mode;
      Component = route.component;
    }

    // wait until the user loads to avoid race conditions
    if (!user) return null;

    // this key allows us to remount the component when the params change
    let key = JSON.stringify(params);

    return (
      <React.Fragment>
        <ModalDialogContainer>
          <VulcanProvider params={params}>
            <ThemeProvider theme={theme}>
              <div id='ui-container'>
                <RemountOnParamsChange params={params}>
                  <Notifications />
                  <VulcanNav environment={environment} mode={mode} />
                  <Messages />
                  <Component key={key} {...params} />
                </RemountOnParamsChange>
              </div>
            </ThemeProvider>
          </VulcanProvider>
        </ModalDialogContainer>
      </React.Fragment>
    );
  }
}

// Used To simulate remount whenever router params changes
// TODO: Can this be removed post refactoring?
function RemountOnParamsChange({params, children}) {
  return children;
}

export default connect(
  (state) => ({location: state.location, user: selectUser(state)}),
  {showMessages, updateLocation}
)(VulcanUI);
