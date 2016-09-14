/**
 *  Barebones XHR/AJAX wrapper
 *  Works like the JQuery AJAX function, but without JQuery
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