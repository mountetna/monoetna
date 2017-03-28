import React from 'react';
import ReactDOM from 'react-dom';
import { Provider } from 'react-redux';

import BasicViewModel from './models/basic-view-model';
import BasicViewContainer from './components/basic-view-container';
import JanusLoggerController from './janus-logger-controller';

class BasicViewController{

  constructor(){

    this['model'] = null
    this['janusLogger'] = new JanusLoggerController();

    this.initDataStore();
    this.buildUI();

    /*
     * We pass in the store since we bind events to it. The AJAX callbacks will
     * be dispatched using the store.
     */
    this['janusLogger']['model']['store'] = this['model']['store'];
    this['janusLogger'].checkLog();
  }

  initDataStore(){

    this['model'] = new BasicViewModel();

    // Event hooks from the UI to the Controller
    this['model']['store'].subscribe(()=>{ 

      var lastAction = this['model']['store'].getState()['lastAction'];
      this.routeAction(lastAction);
    });
  }

  buildUI(){

    ReactDOM.render(

      <Provider store={ this['model']['store'] }>

        <BasicViewContainer />
      </Provider>,
      document.getElementById('ui-group')
    );
  }

  routeAction(action){

    switch(action['type']){

      case 'LOG_IN':

        var email = action['data']['email'];
        var password = action['data']['pass'];
        this['janusLogger'].logIn(email, password);
        break;
      case 'LOG_OUT':

        this['janusLogger'].logOut(COOKIES.getItem(TOKEN_NAME));
        break;
      case 'LOGGED_OUT':
      case 'NOT_LOGGED':

        window.location = LOGGED_OUT_ADDR();
        break;
      default:

        // none
        break;
    }
  }
}

// Initilize the class.
var basicViewController = new BasicViewController();