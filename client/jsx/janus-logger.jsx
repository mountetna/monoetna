/*
 * Contains the functional items that make the Janus client work.
 */
export default class JanusLogger{

  constructor(){

    this['model'] = {};
    this['model']['store'] = null;
  }

  checkLog(callback){

    if(COOKIES.hasItem(TOKEN_NAME)){

      //Serialize the request for POST
      var logItems = 'token='+ COOKIES.getItem(TOKEN_NAME);

      AJAX({

        'url': METIS_ADDR +'/check',
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

    if(response['success'] && response['logged']){

      this.logInResponse(response);
    }
    else if(response['success'] && !response['logged']){

      this.notLogged();
    }
    else{

      this.logError();
    }
  }

  logInResponse(response){

    if(response['success']){

      //Set the token to the cookies so it may be used by multiple UI programs.
      COOKIES.setItem(TOKEN_NAME, response['user_info']['auth_token'], Infinity, '/', 'ucsf.edu');

      //Set the token to the local Redux store.
      var data = { 
        
        'userEmail': response['user_info']['email'],
        'authToken': response['user_info']['auth_token'],
        'firstName': response['user_info']['first_name'],
        'lastName': response['user_info']['last_name'],
        'userId': response['user_info']['user_id'],
        'permissions': response['user_info']['permissions']
      };
      var action = { 'type': 'LOGGED_IN', data: data };
    }
    else{

      var action = { 'type': 'LOG_ERROR' };
    }

    this['model']['store'].dispatch(action);
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

  logOut(callback){

    var state = this['model']['store'].getState();
    var email = state['userInfo']['userEmail'];

    var logItems = [

      'email='+ email,
      'token='+ COOKIES.getItem(TOKEN_NAME)
    ];

    AJAX({

      'url': '/logout',
      'method': 'POST',
      'sendType': 'serial',
      'returnType': 'json',
      'data': logItems.join('&'),
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

    var action = { 'type': 'NOT_LOGGED' };
    this['model']['store'].dispatch(action);
  }

  logError(response){

    var action = { 'type': 'LOG_ERROR' };
    console.log(response);
  }

  ajaxError(xhr, config, error){

    console.log(xhr, config, error);
  }
}