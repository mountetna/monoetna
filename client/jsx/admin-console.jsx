import React from 'react';
import ReactDOM from 'react-dom';
import { Provider } from 'react-redux';

import JanusLogger from './janus-logger';
import KeyboardShortcuts from './keyboard-shortcuts';
import AdminModel from './models/admin-model';
import AdminUIContainer from './components/admin-ui-container';

class AdminConsole{

  constructor(){

    this['model'] = null;
    this['janusLogger'] = new JanusLogger();

    this.initDataStore();

    /*
     * We pass in the store since we bind events to it. The AJAX callbacks will
     * be dispatched using the store.
     */
    this['janusLogger']['model']['store'] = this['model']['store'];
    this['janusLogger'].checkLog();

    this.buildUI();
  }

  buildUI(){

    ReactDOM.render(

      <Provider store={ this['model']['store'] }>

        <AdminUIContainer />
      </Provider>,
      document.getElementById('ui-group')
    );
  }

  initDataStore(){

    this['model'] = new AdminModel();

    // Event hooks from the UI to the Controller
    this['model']['store'].subscribe(()=>{ 

      var lastAction = this['model']['store'].getState()['lastAction'];
      this.routeAction(lastAction);
    });
  }

  routeAction(action){

    switch(action['type']){

      case 'NOT_LOGGED':

        window.location = '/';
        break;

      case 'LOGGED_IN':

        this.fetchAdminData();
        break;
      default:

        // none
        break;
    }
  }

  fetchAdminData(){

    var appState = this['model']['store'].getState()['appState'];
    if(!appState['adminPerms']){

      window.location = '/';
    }
    else{

      this.adminDataCall('/get-projects');
      this.adminDataCall('/get-permissions');
      if(appState['masterPerms']) this.adminDataCall('/get-users');
    }
  }

  adminDataCall(endPoint){

    if(COOKIES.hasItem(TOKEN_NAME)){

      //Serialize the request for POST
      var logItems = 'token='+ COOKIES.getItem(TOKEN_NAME);

      AJAX({

        'url': METIS_ADDR + endPoint,
        'method': 'POST',
        'sendType': 'serial',
        'returnType': 'json',
        'data': logItems,
        'success': this['adminDataResponse'].bind(this),
        'error': this['ajaxError'].bind(this)
      });
    }
  }

  adminDataResponse(response){

    console.log(response);
    var action = null;
    if(response['success']){

      if(response.hasOwnProperty('projects')){

        action = { 'type': 'HAS_PROJECTS', 'data': response['projects'] };
      }
      else if(response.hasOwnProperty('permissions')){

        action = { 'type': 'HAS_PERMISSIONS', 'data': response['permissions'] };
      }
      else if(response.hasOwnProperty('users')){

        action = { 'type': 'HAS_USERS', 'data': response['users'] };
      }
      else{

        action = null;
      }
    }

    if(action != null) this['model']['store'].dispatch(action);
  }

  ajaxError(xhr, config, error){

    console.log(xhr, config, error);
  }
}

// Initilize the class.
var adminConsole = new AdminConsole();