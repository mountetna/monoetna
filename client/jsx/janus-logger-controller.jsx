/*
 * Contains the functional items that make the Janus client work.
 */
export default class JanusLoggerController{

  constructor(){

    this['model'] = {};
    this['model']['store'] = null;
  }

  checkLog(callback){

    if(COOKIES.hasItem(TOKEN_NAME)){

      //Serialize the request for POST
      var logItems = 'token='+ COOKIES.getItem(TOKEN_NAME);

      AJAX({

        'url': '/check',
        'method': 'POST',
        'sendType': 'serial',
        'returnType': 'json',
        'data': logItems,
        'success': this['checkLogResponse'].bind(this),
        'error': this['ajaxError'].bind(this)
      });
    }
    else{

      this.notLogged();
    }
  }

  checkLogResponse(response){

    if(response['success'] && response['user_info']){

      this.logInResponse(response);
    }
    else{

      this.notLogged();
    }
  }

  logIn(email, pass, callback){

    var logItems = [

      'email='+ email,
      'pass='+ pass
    ];

    AJAX({

      'url': '/login',
      'method': 'POST',
      'sendType': 'serial',
      'returnType': 'json',
      'data': logItems.join('&'),
      'success': this['logInResponse'].bind(this),
      'error': this['ajaxError'].bind(this)
    });
  }

  logInResponse(response){

    if(response['success']){

      // Set the token to the cookies so it may be used by multiple UI programs.
      var token = response['user_info']['token'];
      COOKIES.setItem(TOKEN_NAME, token, Infinity, '/', 'ucsf.edu');

      // Check to see if the user has any permissions.
      if(response['user_info']['permissions']['length'] == 0){

        this.setNoPermissionAlert();
        return;
      }

      var action = {

        'type': 'LOGGED_IN',
        'data': this.formLoginResponseData(response)
      };
    }
    else{

      var action = { 'type': 'LOG_ERROR' };
    }

    this['model']['store'].dispatch(action);
  }

  setNoPermissionAlert(){

    var msg = 'It looks like you have no permissions in our system, '
    msg += 'which means you cannot permform any tasks in our system. '
    msg += 'Please contact the system administrator to get some permissions.';

    alert(msg);

    var token = response['user_info']['auth_token'];
    this.logOut(token);
  }

  formLoginResponseData(response){

    return {

      'userEmail': response['user_info']['email'],
      'authToken': response['user_info']['token'],
      'firstName': response['user_info']['first_name'],
      'lastName': response['user_info']['last_name'],
      'userId': response['user_info']['user_id'],
      'permissions': response['user_info']['permissions']
    };
  }

  logOut(token){

    AJAX({

      'url': '/logout',
      'method': 'POST',
      'sendType': 'serial',
      'returnType': 'json',
      'data': 'token='+ token,
      'success': this['logOutResponse'].bind(this),
      'error': this['ajaxError'].bind(this)
    });
  }

  logOutResponse(response){

    if(response['success'] && !response['logged']){

      COOKIES.removeItem(TOKEN_NAME, '/', 'ucsf.edu');
      var action = { 'type': 'LOGGED_OUT' };
      this['model']['store'].dispatch(action);
    }
    else{

      this.logError();
    }
  }

  notLogged(){

    if(COOKIES.hasItem(TOKEN_NAME)){

      COOKIES.removeItem(TOKEN_NAME, '/', 'ucsf.edu');
    }

    var action = { 'type': 'NOT_LOGGED' };
    this['model']['store'].dispatch(action);
  }

  logError(response){

    if(COOKIES.hasItem(TOKEN_NAME)){

      COOKIES.removeItem(TOKEN_NAME, '/', 'ucsf.edu');
    }

    var action = { 'type': 'LOG_ERROR' };
    this['model']['store'].dispatch(action);
  }

  ajaxError(xhr, config, error){

    console.log(xhr, config, error);
  }
}