import * as Redux from 'redux';

import AdminReducer from './admin-reducer';
import JanusLogReducer from './janus-log-reducer';
import LastActionReducer from './last-action-reducer';

export default class AdminModel{

  constructor(){

    var adminReducer = new AdminReducer();
    var janusLogReducer = new JanusLogReducer();
    var lastAction = new LastActionReducer();
    var reducer = Redux.combineReducers({

      'adminInfo': adminReducer.reducer(),
      'userInfo': janusLogReducer.reducer(),
      'lastAction': lastAction.reducer()
    });

    var defaultState = {

      'adminInfo': {

        'users': [],
        'projects': [],
        'permissions': []
      },

      'userInfo': {

        'userEmail': '',
        'authToken': '',
        'firstName': '',
        'lastName': '',
        'permissions': [],

        'masterPerms': false,

        'loginStatus': false,
        'loginError': false,
        'loginErrorMsg': 'Invalid sign in.'
      }
    };

    this.store = Redux.createStore(reducer, defaultState);
  }
}