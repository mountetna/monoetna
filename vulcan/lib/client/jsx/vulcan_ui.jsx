import * as React from 'react';
import {connect} from 'react-redux';
import {findRoute, setRoutes} from './router';

// Components.
import Browser from './components/workflow/browser';
import RootView from 'etna-js/components/RootView';
import VulcanNav from './components/vulcan_nav';
import Messages from 'etna-js/components/messages';
import {selectUser} from 'etna-js/selectors/user-selector';

import {showMessages} from 'etna-js/actions/message_actions';
import {updateLocation} from 'etna-js/actions/location_actions';

import {ModalDialogContainer} from 'etna-js/components/ModalDialogContainer';
import {Notifications} from 'etna-js/components/Notifications';

const ROUTES = [
  {
    template: '',
    component: RootView,
    mode: 'home'
  },
  {
    template: ':project_name/',
    component: Browser,
    mode: 'workflow'
  },
  {
    name: 'workflow',
    template: ':project_name/workflow/',
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
