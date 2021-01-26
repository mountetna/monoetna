import * as React from 'react';
import {connect} from 'react-redux';
import {findRoute, setRoutes} from './router';

// Components.
import WorkflowBrowser from './components/workflow/browser';
import RootView from 'etna-js/components/RootView';
import VulcanNav from './components/vulcan_nav';
import Messages from './components/messages';

import {showMessages} from './actions/message_actions';
import {updateLocation} from './actions/location_actions';

import {selectUser} from './selectors/user_selector';
import {ModalDialogContainer} from 'etna-js/components/ModalDialogContainer';
import {Notifications} from 'etna-js/components/Notifications';

const ROUTES = [
  {
    template: '',
    component: RootView,
    mode: 'home'
  },
  {
    name: 'workflow',
    template: 'workflow',
    component: WorkflowBrowser,
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

  render() {
    let {location, showMessages, environment, user} = this.props;

    let {route, params} = findRoute(location, ROUTES);
    let Component;
    let mode;

    if (!route) {
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
          <div id='ui-container'>
            <Notifications />
            <VulcanNav environment={environment} mode={mode} />
            <Messages />
            <Component key={key} {...params} />
          </div>
        </ModalDialogContainer>
      </React.Fragment>
    );
  }
}

export default connect(
  (state) => ({location: state.location, user: selectUser(state)}),
  {showMessages, updateLocation}
)(VulcanUI);
