import * as Redux from 'redux';

import AdminReducer from './admin-reducer';
import LastActionReducer from './last-action-reducer';

export default class AdminModel{

  constructor(){

    var appState = new AdminReducer();
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

        'adminPerms': false,
        'masterPerms': false,

        'loginStatus': false,
        'loginError': false,
        'loginErrorMsg': 'Invalid sign in.',

        'users': [],
        'projects': [],
        'permissions': []
      }
    };

    this.store = Redux.createStore(reducer, defaultState);
  }
}