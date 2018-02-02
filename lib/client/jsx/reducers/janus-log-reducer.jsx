export default class JanusLogReducer{
  reducer(){
    return (state={}, action)=>{
      switch(action.type){
        case 'LOG_IN':
          return {
            ...state,
            userEmail: action.data.email
          }

        case 'LOGGED_IN':
          let userInfo = Object.assign({}, state);

          /* 
           * Copy the new data from the auth server to the local Redux store.
           * Also keep an eye on that 'cleanPermissions' function. The client
           * wants it's vars in camel case.
           */
          for(let key in action.data){
            let userItem = action.data[key];
            if(key == 'permissions') userItem = this.cleanPermissions(userItem);
            userInfo[key] = userItem;
          }

          userInfo.loginStatus = true;
          userInfo.loginError = false;
          let perms = userInfo.permissions;
          userInfo.masterPerms = this.checkAdminPermissions(perms);
          return userInfo;

        case 'LOGGED_OUT':
          return {
            userId: null,
            permissions: [],
            masterPerms: false,
            loginStatus: false,
            logError: false,
            loginErrorMsg: 'Invalid sign in.'
          };

        case 'LOG_ERROR':
          return {
            ...state,
            loginStatus: false,
            loginError: true
          };

        default:
          return state;
      }
    };
  }

  cleanPermissions(perms){
    for(let index in perms){
      for(let key in perms[index]){
        perms[index][CAMEL_CASE_IT(key)] = perms[index][key];
        if(key.indexOf('_') != -1) delete perms[index][key];
      }
    }

    return perms;
  }

  checkAdminPermissions(perms){
    // Check for administration privileges.
    let masterPerms = false;
    for(let index in perms){
      if(perms[index].role == 'administrator'){
        if(perms[index].projectName == 'administration'){
          masterPerms = true;
        }
      }
    }

    return masterPerms;
  }
}
