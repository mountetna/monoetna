import * as Redux from 'redux';

import MetisReducer from './metis-reducer';
import LastActionReducer from './last-action-reducer';

export default class MetisModel{

  constructor(){

    var appState = new MetisReducer();
    var lastAction = new LastActionReducer();
    var reducer = Redux.combineReducers({

      'appState': appState.reducer(),
      'lastAction': lastAction.reducer()
    });

    var defaultState = {

      'appState': {

        'userInfo': {

          'userEmail': '',
          'authToken': '',
          'firstName': '',
          'lastName': '',
          'permissions': []
        },

        'fileList': [],
        'fileUploads': [],

        'loginStatus': false,
        'loginError': false,
        'loginErrorMsg': 'Invalid sign in.'
      }
    };

    this.store = Redux.createStore(reducer, defaultState);
  }
}