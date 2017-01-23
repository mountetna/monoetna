import * as Redux from 'redux';

import MetisReducer from './metis-reducer';
import JanusLogReducer from './janus-log-reducer';
import LastActionReducer from './last-action-reducer';

export default class MetisModel{

  constructor(){

    var metisReducer = new MetisReducer();
    var janusLogReducer = new JanusLogReducer();
    var lastAction = new LastActionReducer();
    var reducer = Redux.combineReducers({

      'fileData': metisReducer.reducer(),
      'userInfo': janusLogReducer.reducer(),
      'lastAction': lastAction.reducer()
    });

    var defaultState = {

      'fileData': {

        'fileList': [],
        'fileUploads': [],
        'fileFails': []
      },

      'userInfo': {

        'userId': null,
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