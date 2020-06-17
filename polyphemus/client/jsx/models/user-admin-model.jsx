import * as Redux from 'redux';

import UserAdminReducer from './user-admin-reducer';
import JanusLogReducer from './janus-log-reducer';
import LastActionReducer from './last-action-reducer';

export default class UserAdminInfo{

  constructor(){

    var userAdminReducer = new UserAdminReducer();
    var janusLogReducer = new JanusLogReducer();
    var lastAction = new LastActionReducer();
    var reducer = Redux.combineReducers({

      'adminInfo': userAdminReducer.reducer(),
      'userInfo': janusLogReducer.reducer(),
      'lastAction': lastAction.reducer()
    });

    var defaultState = {

      'adminInfo': {

        'users': [],
        'projects': [],
        'permissions': [] // These are all the permissions in the system.
      },

      'userInfo': {

        'userEmail': '',
        'authToken': '',
        'firstName': '',
        'lastName': '',
        'permissions': [], // These are the permissions for the logged user.

        'masterPerms': false,

        'loginStatus': false,
        'loginError': false,
        'loginErrorMsg': 'Invalid sign in.'
      }
    };

    this.store = Redux.createStore(reducer, defaultState);
  }
}