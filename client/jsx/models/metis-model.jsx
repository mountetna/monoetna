import * as Redux from 'redux';

import MetisReducer from './metis-reducer';
import LastActionReducer from './last-action-reducer';

export default class MetisModel{

  constructor(){

    var metisState = new MetisReducer();
    var lastAction = new LastActionReducer();
    var reducer = Redux.combineReducers({

      metisState: metisState.reducer(),
      lastAction: lastAction.reducer()
    });

    var defaultState = {

      metisState: {

        fileList: [],
        fileUploads: [],
        userInfo: {

          user_email: 'jasondcater@gmail.com',
          authorization_token: 'blahf8cfc63531b6ed753bd536f6e12d578c'
        }
      }
    };

    this.store = Redux.createStore(reducer, defaultState);
  }
}