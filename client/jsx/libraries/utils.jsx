/*
 * Global Variables
 */
var BLOB_SIZE = 1024; // in bytes
var NEXT_BLOB_SIZE = 1024; // in bytes
var MIN_BLOB_SIZE = 1024; // in bytes

/*
 * In milliseconds, the amount of time to transfer one blob. This ultimately
 * sets the blob size.
 */
var TRANSFER_TIME = 2000; 

var TOKEN_NAME = 'UCSF_ETNA_AUTH_TOKEN';
var METIS_ADDR = 'http://metis-dev.ucsf.edu';

/*
 * These are the items that are passed between the server and client.
 */
var STATUS_ITEMS = {

  'directory': String,
  'expires': Number,
  'signingAlgorithm': String,
  'hashingAlgorithm': String,
  'startTimestamp': Number,
  'finishTimestamp': Number,
  'authToken': String,
  'authorizationToken': String,
  'token': String,
  'originalName': String,
  'fileName': String,
  'fileSize': Number,
  'userEmail': String,
  'userId': Number,
  'projectId': Number,
  'projectName': String,
  'groupId': Number,
  'groupName': String,
  'dbIndex': String,
  'status': String,
  'role': String,

  'hmacSignature': String,
  
  'currentBlobSize': Number,
  'currentBytePosition': Number,
  'nextBlobHash': String,
  'nextBlobSize': Number,

  'uploadSpeed': Number
};

/* 
 * Generates a psudo random key for React components and data element sorting.
 * DO NOT USE FOR SECURITY!
 */
var GENERATE_RAND_KEY = function(){

  var randKey = '';
  for(var a = 0; a < 4; ++a){

    randKey += (Math.random() * 0xFFFFFF << 0).toString(16);
  }

  return randKey;
}

/*
 * Change a timestamp in to a human readable format.
 */
var PARSE_TIMESTAMP = function(timestamp){

  timestamp = parseInt(timestamp);
  if(isNaN(timestamp) || timestamp < 0){

    timestamp = new Date().getTime() / 1000;
  }

  var dt = new Date(timestamp * 1000);

  var hr = dt.getHours() + 1;
  var mn = dt.getMinutes();

  var mo = dt.getMonth() + 1;
  var dy = dt.getDate();
  var yr = dt.getYear();

  return hr +':'+ mn +', '+ mo +'/'+ dy;
}

/*
 * Change an integer of bytes into a human readable format.  
 * http://stackoverflow.com/questions/10420352/converting-file-size-in-bytes-to-human-readable
 */
var PARSE_BYTES = function(bytes, si){

  var thresh = si ? 1000 : 1024;
  if(Math.abs(bytes) < thresh){

      return bytes + ' B';
  }

  var units = si
      ? ['kB','MB','GB','TB','PB','EB','ZB','YB']
      : ['KiB','MiB','GiB','TiB','PiB','EiB','ZiB','YiB'];

  var u = -1;
  do{

    bytes /= thresh;
    ++u;
  }
  while(Math.abs(bytes) >= thresh && u < units.length - 1);

  return bytes.toFixed(1)+' '+units[u];
}

/*
 * Data from the controller to the worker cannot be passed through messaging
 * while attached to the File Object (the file request/status data is stored on 
 * the file object in redux). We need to externalize the request/status data
 * from the File Object and pass it over 'worker messaging' as it's own 
 * object/hash.
 */
var PARSE_REQUEST = function(file){

  var request = {};

  for(var key in STATUS_ITEMS){

    if(key in file){

      switch(STATUS_ITEMS[key]){

        case Number:

          request[key] = parseInt(file[key]);
          break;
        case String:

          request[key] = String(file[key]);
          break;
        default:

          request[key] = file[key];
          break;
      }
    }
  }

  return request;
}

/*
 * Verification of types and structures of the server responses. This will also
 * transform any strings that are supposed to be integers. 
 * 
 * This function needs to be simplified.
 */

var NORMILIZE_RESPONSE = function(response){

  if('success' in response){

    var success = response['success'];
    if((success == 'true') || (success == 'True') || (success == true)){

      response['success'] == true;
    }
    else if((success == 'false') || (success == 'False') || (success == false)){

      response['success'] == false;
    }
    else{

      return false;
    }
  }

  return response;
}

var TRANSFORM_RESPONSE = function(response){

  var error = false;

  if('request' in response){

    var type = Object.prototype.toString.call(response['request']);
    if(type == '[object String]'){

      response['request'] = JSON.parse(response['request']);
    }
  }

  for(var key in response['request']){

    if(key in STATUS_ITEMS){

      switch(STATUS_ITEMS[key]){

        case String:

          response['request'][key] = String(response['request'][key]);
          break;
        case Number:
          
          response['request'][key] = parseInt(response['request'][key]);
          if(isNaN(response['request'][key])) error = true;
          break;
        default:

          //none
          break;
      }
    }
  }

  return (error) ? error : response;
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

/*
 * Dump a blob onto the console.
 */
var HEX_DUMP = function(blob, fromByte){

  var topLine = ['       00 01 02 03 04 05 06 07 08 09 0A 0B 0C 0D 0E 0F'];
  var underline = '      ';
  for(var a = 0; a < 48; ++a){

    underline += '_';
  }

  var fileReader = new FileReader();

  fileReader.onload = function(progressEvent){

    var arrayBufferNew = this.result;
    var uint8Array = new Uint8Array(arrayBufferNew);

    var hexCount = 0;
    var hexLine = '';
    var hexRows = [];

    for(var b = 0; b < (fromByte%16); ++b){

      hexLine += '   ';
      ++hexCount;
    }

    for(var b = 0; b < uint8Array.length; ++b){

      var hexChar = uint8Array[b].toString(16);
      if(hexChar.length === 1) hexChar = '0'+ hexChar;
      hexLine += hexChar +' ';
      ++hexCount;

      if(hexCount%16 === 0){

        hexRows.push(hexLine);
        hexLine = '';
        hexCount = 0;
      }
    }

    console.log();
    console.log(topLine.join(' '));
    console.log(underline);

    var rowOffset = fromByte - (fromByte%16);

    for(var row in hexRows){

      var offset = (row*16)+rowOffset;
      var lineNum = ('000' + offset.toString(16)).slice(-4);
      console.log(lineNum +' | '+ hexRows[row]);
    }
  }

  fileReader.readAsArrayBuffer(blob);
}

/*  
 * From a file called cookie.js
 *
 * A complete cookies reader/writer framework with full unicode support.
 *
 * https://developer.mozilla.org/en-US/docs/DOM/document.cookie
 *
 * This framework is released under the GNU Public License, version 3 or later.
 * http://www.gnu.org/licenses/gpl-3.0-standalone.html
 *
 * Syntaxes:
 *
 * Cookies.setItem(name, value[, end[, path[, domain[, secure]]]])
 * Cookies.getItem(name)
 * Cookies.removeItem(name[, path], domain)
 * Cookies.hasItem(name)
 * Cookies.keys()
 *
 */

var COOKIES = {

  getItem: function (sKey){

    return decodeURIComponent(document.cookie.replace(new RegExp("(?:(?:^|.*;)\\s*" + encodeURIComponent(sKey).replace(/[\-\.\+\*]/g, "\\$&") + "\\s*\\=\\s*([^;]*).*$)|^.*$"), "$1")) || null;
  },
  setItem: function(sKey, sValue, vEnd, sPath, sDomain, bSecure){

    if(!sKey || /^(?:expires|max\-age|path|domain|secure)$/i.test(sKey)){ return false; }
    
    var sExpires = "";
    if(vEnd){

      switch(vEnd.constructor){

        case Number:

          sExpires = vEnd === Infinity ? "; expires=Fri, 31 Dec 9999 23:59:59 GMT" : "; max-age=" + vEnd;
          break;
        case String:

          sExpires = "; expires=" + vEnd;
          break;
        case Date:

          sExpires = "; expires=" + vEnd.toUTCString();
          break;
      }
    }

    document.cookie = encodeURIComponent(sKey) + "=" + encodeURIComponent(sValue) + sExpires + (sDomain ? "; domain=" + sDomain : "") + (sPath ? "; path=" + sPath : "") + (bSecure ? "; secure" : "");
    return true;
  },
  removeItem: function(sKey, sPath, sDomain){

    if(!sKey || !this.hasItem(sKey)){ return false; }

    document.cookie = encodeURIComponent(sKey) + "=; expires=Thu, 01 Jan 1970 00:00:00 GMT" + ( sDomain ? "; domain=" + sDomain : "") + ( sPath ? "; path=" + sPath : "");
    return true;
  },
  hasItem: function(sKey){

    return (new RegExp("(?:^|;\\s*)" + encodeURIComponent(sKey).replace(/[\-\.\+\*]/g, "\\$&") + "\\s*\\=")).test(document.cookie);
  },

  /* optional method: you can safely remove it! */ 
  keys: function(){

    var aKeys = document.cookie.replace(/((?:^|\s*;)[^\=]+)(?=;|$)|^\s*|\s*(?:\=[^;]*)?(?:\1|$)/g, "").split(/\s*(?:\=[^;]*)?;\s*/);
    for(var nIdx = 0; nIdx < aKeys.length; nIdx++) { aKeys[nIdx] = decodeURIComponent(aKeys[nIdx]); }

    if(aKeys.length == 1){

      if(aKeys[0] == ""){

        return null;
      }
    }
    
    return aKeys;
  }
};

/*
 * This is only to prevent sending excessive queries to the server. There is
 * server side email validation too.
 */
var VALIDATE_EMAIL = function(email){

    var re = /^(([^<>()\[\]\\.,;:\s@"]+(\.[^<>()\[\]\\.,;:\s@"]+)*)|(".+"))@((\[[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}])|(([a-zA-Z\-0-9]+\.)+[a-zA-Z]{2,}))$/;
    return re.test(email);
}

/*
 * In Ruby (the server side) variables/hash keys are named in snake_case.
 * Here in JS (the client side) we are using camelCase. 
 * This is a little util to transform snake_case to camelCase and back again.
 */
var CAMEL_CASE_IT = function(string){

  var regexp = /_/g;
  var match, matches = [];
  while((match = regexp.exec(string)) != null){

    matches.push(match.index);
  }

  // We should start from the back of the string and replace going forward.
  matches.reverse();
  for(var index in matches){

    var ind = matches[index];
    var fromCharacter = '_'+string.charAt(ind+1);
    var toCharacter = string.charAt(ind+1).toUpperCase()
    string = string.replace(fromCharacter, toCharacter);
  }

  return string;
}

var SNAKE_CASE_IT = function(str){

  for(var a = 0; a < str['length']; ++a){

    var chr = str[a];
    if(chr != chr.toLowerCase() && a != 0 && a != (str['length']-1)){

      str = str.substr(0,a) +'_'+ chr.toLowerCase() + str.substr(a+1);
    }
  }

  return str;
}

// After we log out of Janus we must log out of Shibboleth.
var LOGGED_OUT_ADDR = function(){

  var base = 'https://janus-stage.ucsf.edu/Shibboleth.sso/Logout';
  return base+'?return=http%3A%2F%2Fmetis-dev.ucsf.edu%2Flogged-out';
}

var NOT_LOGGED_ADDR = function(){

  var base = 'https://janus-stage.ucsf.edu/login';
  return base+'?refer=http%3A%2F%2Fmetis-dev.ucsf.edu%2F';
}