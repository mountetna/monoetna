/*
 * Global Variables
 */
var BLOB_SIZE = 1024000; // in bytes

/*
 * Verification of types and structures of the server responses. This will also
 * transform any strings that are supposed to be integers. 
 * 
 * This function needs to be simplified.
 */
var VERIFY_AND_TRANSFORM = function(response){

  if('request' in response){

    response['request'] = JSON.parse(response['request']);
  }
  else{

    return false;
  }

  if('byte_count' in response){
    
    response['byte_count'] = parseInt(response['byte_count']);
  }
  else{

    return false;
  }

  if('signature' in response){

    response['signature'] = String(response['signature']);
  }
  else{

    return false;
  }

  if('status' in response){

    response['status'] = String(response['status']);
  }
  else{

    return false;
  }

  if('success' in response){

    if((response['success'] == 'true') || 
       (response['success'] == 'True') ||
       (response['success'] == true)){

      response['success'] == true;
    }
    else if((response['success'] == 'false') || 
       (response['success'] == 'False') ||
       (response['success'] == false)){

      response['success'] == false;
    }
    else{

      return false
    }
  }

  var error = false;
  var returnItems = [

    {algorithm: String},
    {authorization_token: String},
    {current_blob_size: Number},
    {current_byte_position: Number},
    {directory: String},
    {expires: Number},
    {file_name: String},
    {file_size: Number},
    {next_blob_hash: String},
    {next_blob_size: Number},
    {original_name: String},
    {signature: String},
    {timestamp: Number},
    {user_email: String}
  ]

  for(var index in returnItems){

    for(var key in returnItems[index]){

      var value = returnItems[index][key];
      if(key in response['request']){

        switch(value){

          case String:

            response['request'][key] = String(response['request'][key]);
            break;
          case Number:
            
            response['request'][key] = parseInt(response['request'][key]);
            break;
          default:

            //none
            break;
        }
      }
      else{

        error = true;
      }
    }
  }

  if(error){

    return false;
  }
  else{

    return response;
  }
}

/*
 * Barebones XHR/AJAX wrapper
 * Works like the JQuery AJAX function, but without JQuery
 */
var AJAX = function(config){

  var url = config.url;
  var method = config.method;
  var sendType = config.sendType;
  var returnType = config.returnType;
  var success = config.success;
  var error = config.error;
  var data = (config.data === undefined) ? "" : config.data;

  var xhr = new XMLHttpRequest();
  
  /*
   * Response Section
   */
  xhr.onreadystatechange = function(){

    switch(xhr.readyState){

      case 0:

        //none
        break;
      case 1:

        //connection opened
        break;
      case 2:

        //headers received
        var type = xhr.getResponseHeader('Content-Type');
        break;
      case 3:

        //loading
        break;
      case 4:

        if(xhr.status === 200){

          switch(returnType.toLowerCase()){

            case 'json':

              success(JSON.parse(xhr.responseText));
              break;
            case 'html':

              success(elem);
              break;  
            default:

              success(xhr.responseText);
              break;
          }
        }
        else{

          error(xhr, config, 'error');
        }
        break;
      default:
        break;
    }
  };
  
  /*
   * Execution/Call Section
   */
  xhr.open(method, url, true);

  switch(method.toLowerCase()){

    case 'post':

      /*
       * If we are using a FormData object then the header is already 
       * set appropriately
       */
      if(sendType.toLowerCase() != 'file'){

        var headerType = 'application/x-www-form-urlencoded';
        xhr.setRequestHeader('Content-type', headerType);
      }
      xhr.send(data);
      break;
    case 'get':

      xhr.send();
      break;
    default:

      error(xhr, config, 'Unknown HTTP Method : "'+ method.toLowerCase() +'"');
      break;
  }
}