export default class JanusLogReducer{

  reducer(){

    return (state={}, action)=>{

      switch(action['type']){

        case 'LOG_IN':

          var userInfo = Object.assign({}, state);
          userInfo['userEmail'] = action['data']['email'];
          return userInfo;

        case 'LOGGED_IN':

          var userInfo = Object.assign({}, state);

          /* 
           * Copy the new data from the auth server to the local Redux store.
           * Also keep an eye on that 'cleanPermissions' function. The client
           * wants it's vars in camel case.
           */
          for(var key in action['data']){

            var userItem = action['data'][key];
            if(key == 'permissions') userItem = this.cleanPermissions(userItem);
            userInfo[key] = userItem;
          }

          userInfo['loginStatus'] = true;
          userInfo['loginError'] = false;
          var perms = userInfo['permissions'];
          userInfo['masterPerms'] = this.checkAdminPermissions(perms);
          return userInfo;

        case 'LOGGED_OUT':

          var userInfo = Object.assign({}, state);

          // Clear the local data.
          for(var key in userInfo){

            userInfo[key] = '';
          }

          userInfo['permissions'] = [];
          userInfo['masterPerms'] = false;
          userInfo['loginStatus'] = false;
          userInfo['logError'] = false;
          userInfo['loginErrorMsg'] = 'Invalid sign in.';
          return userInfo;
        case 'LOG_ERROR':

          var userInfo = Object.assign({}, state);
          userInfo['loginStatus'] = false;
          userInfo['loginError'] = true;
          return userInfo;

        default:

          var userInfo = Object.assign({}, state);
          return userInfo;
      }
    };
  }

  cleanPermissions(perms){

    for(var index in perms){

      for(var key in perms[index]){

        perms[index][CAMEL_CASE_IT(key)] = perms[index][key];
        if(key.indexOf('_') != -1) delete perms[index][key];
      }
    }

    return perms;
  }

  checkAdminPermissions(perms){

    // Check for administration privileges.
    var masterPerms = false;
    for(var index in perms){

      if(perms[index]['role'] == 'administrator'){

        if(perms[index]['projectName'] == 'administration'){

          masterPerms = true;
        }
      }
    }

    return masterPerms;
  }
}