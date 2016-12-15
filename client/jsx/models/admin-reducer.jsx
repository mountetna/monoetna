export default class AdminReducer{

  reducer(){

    return (state = {}, action)=>{

      switch(action['type']){

        case 'LOG_IN':

          var nextState = Object.assign({}, state);
          nextState['userInfo']['userEmail'] = action['data']['email'];
          return nextState;

        case 'LOGGED_IN':

          var nextState = Object.assign({}, state);

          /* Copy the new data from the auth server to the local Redux store.
           * Also keep an eye on that 'cleanPermissions' function. The client
           * wants it's vars in camel case.
           */
          for(var key in action['data']){

            var userItem = action['data'][key];
            if(key == 'permissions') userItem = this.camelCaseIt(userItem);
            nextState['userInfo'][key] = userItem;
          }

          // Set the login state.
          nextState['loginStatus'] = true;
          nextState['loginError'] = false;

          // Set the admin permssions.
          var perms = nextState['userInfo']['permissions'];
          var adminRights = this.checkAdminPermissions(perms);
          nextState['adminPerms'] = adminRights[0];
          nextState['masterPerms'] = adminRights[1];

          return nextState;
        case 'LOGGED_OUT':

          var nextState = Object.assign({}, state);

          for(var key in nextState['userInfo']){

            nextState['userInfo'][key] = '';
          }

          nextState['loginStatus'] = false;
          nextState['logError'] = false;

          return nextState;
        case 'LOG_ERROR':

          var nextState = Object.assign({}, state);
          nextState['loginStatus'] = false;
          nextState['loginError'] = true;
          nextState['loginErrorMsg'] = 'Invalid sign in.';
          return nextState;

        case 'HAS_PROJECTS':

          var nextState = Object.assign({}, state);
          nextState['projects'] = this.camelCaseIt(action['data']);
          return nextState;

        case 'HAS_PERMISSIONS':

          var nextState = Object.assign({}, state);
          nextState['permissions'] = this.camelCaseIt(action['data']);

          // Generate some random id's for react keys.
          for(var a = 0; a < nextState['permissions']['length']; ++a){

            nextState['permissions'][a]['reactKey'] = GENERATE_RAND_KEY()
          }

          return nextState;

        case 'ADD_PERMISSION':

          var nextState = Object.assign({}, state);

          var randKey = GENERATE_RAND_KEY()
          var newPermission = {
  
            'id': 'permission-'+randKey,
            'projectId': null,
            'projectName': '',
            'role': '',
            'userEmail': '',
            'userId': null,
            'reactKey': randKey
          }

          nextState['permissions'].unshift(newPermission);
          return nextState;

        case 'REMOVE_UNSAVED_PERMISSION':

          var nextState = Object.assign({}, state);

          for(var a = 0; a < nextState['permissions']['length']; ++a){

            if(nextState['permissions'][a]['id'] == action['permissionId']){

              nextState['permissions'].splice(a, 1);
            }
          }
          return nextState;

        case 'HAS_USERS':

          var nextState = Object.assign({}, state);
          nextState['users'] = this.camelCaseIt(action['data']);
          return nextState;
          
        default:

          var nextState = Object.assign({}, state);
          return nextState;
      }
    };
  }

  camelCaseIt(object){

    for(var index in object){

      for(var key in object[index]){

        object[index][CAMEL_CASE_IT(key)] = object[index][key];
        if(key.indexOf('_') != -1) delete object[index][key];
      }
    }

    return object;
  }

  checkAdminPermissions(perms){

    // Check for administration privileges.
    var adminPerms = false;
    var masterPerms = false;
    for(var index in perms){

      if(perms[index]['role'] == 'administrator'){

        adminPerms = true;

        if(perms[index]['projectName'] == 'administration'){

          masterPerms = true;
        }
      }
    }

    return [adminPerms, masterPerms];
  }
}