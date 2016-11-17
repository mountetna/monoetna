import * as Redux from 'redux';

import MetisReducer from './metis-reducer';
import LastActionReducer from './last-action-reducer';

export default class MetisModel{

  constructor(){

    var appState = new MetisReducer();
    var lastAction = new LastActionReducer();
    var reducer = Redux.combineReducers({

      appState: appState.reducer(),
      lastAction: lastAction.reducer()
    });

    var defaultState = {

      appState: {

        userInfo: {

          user_email: '',
          authorization_token: '',
          first_name: '',
          last_name: '',

          userEmail: '',
          authToken: '',
          firstName: '',
          lastName: ''
        },

        fileList: [],
        fileUploads: [],

        loginStatus: false,
        loginError: false
      }
    };

    this.store = Redux.createStore(reducer, defaultState);
  }
}