export default class MetisReducer{

  reducer(){

    return (state = {}, action)=>{

      switch(action['type']){

        case 'FILE_SELECTED':
          
          var fileData = Object.assign({}, state);

          // MOD START
          var fileObject = action['fileObject'];
          fileObject['fileName'] = fileObject['name'];
          fileObject['originalName'] = fileObject['name'];
          fileObject['fileSize'] = fileObject['size'];
          fileObject['currentBytePosition'] = 0;
          fileObject['status'] = 'unauthorized';

          /*
           * This item will be replaced by a Redis id. We will use it in Redis
           * as a unique entry identifier and also in React as a DOM key.
           * For new uploads we just generate a psudo random string, here at the
           * client, to hold a place.
           */
          fileObject['redisIndex'] = GENERATE_RAND_KEY(); 
          // MOD END
          
          fileData['fileUploads'].push(fileObject);
          return fileData;

        case 'FILE_UPLOAD_AUTHORIZED':

          var fileData = Object.assign({}, state);
          var fileUploads = fileData['fileUploads'];

          var authResponse = Object.assign({}, action['authResponse']);
          var fileUpload = null;
          var fileUploadIndex = 0;

          // Select the file to upload from the redux store using the old index.
          var oldIndex = authResponse['request']['old_index'];
          for(var a = 0; a < fileUploads.length; ++a){

            if(fileUploads[a]['redisIndex'] == oldIndex){

              fileUploadIndex = a; 
              fileUpload = fileUploads[a];
              break;
            }
          }
          delete authResponse['request']['old_index'];

          // Append the signature and set the server current byte to 0
          fileUpload['signature'] = authResponse['signature'];
          fileUpload['currentBytePosition'] = 0;

          // Append all of the request items to the local file object.
          authResponse['request'] = this.camelCaseIt(authResponse['request']);
          fileUpload = Object.assign(fileUpload, authResponse['request']);

          fileData['fileUploads'][fileUploadIndex] = fileUpload;
          return fileData;

        case 'FILE_INITIALIZED':

          var fileData = Object.assign({}, state);
          var fileUploads = fileData['fileUploads'];
          var initResponse = Object.assign({}, action['initResponse']);
          initResponse['request'] = this.camelCaseIt(initResponse['request']);
          var request = initResponse['request'];

          for(var a = 0; a < fileUploads['length']; ++a){

            if(fileUploads[a]['redisIndex'] == request['redisIndex']){

              //uploadFile = fileUploads[a];
              for(var key in request){

                fileUploads[a][key] = request[key];
              }
              break;
            }
          }
          return fileData;

        case 'FILE_UPLOAD_ACTIVE':

          var fileData = Object.assign({}, state);
          var fileUploads = fileData['fileUploads'];

          var response = action['uploadResponse'];
          var index = response['request']['redis_index'];
          var fileUpload = null;
          var fileUploadIndex = 0;

          // Select the file to upload from the redux store.
          for(var a = 0; a < fileUploads.length; ++a){

            if(fileUploads[a]['redisIndex'] == index){

              fileUploadIndex = a; 
              fileUpload = fileUploads[a];
              break;
            }
          }

          // Append all of the request items to the local file object.
          response['request'] = this.camelCaseIt(response['request']);
          fileUpload = Object.assign(fileUpload, response['request']);

          fileData['fileUploads'][fileUploadIndex] = fileUpload;
          return fileData;

        case 'FILE_UPLOAD_COMPLETE':

          var fileData = Object.assign({}, state);
          var fileUploads = fileData['fileUploads'];
          var fileList = fileData['fileList'];

          // MOD START
          var result = this.camelCaseIt(action['uploadResponse']['result']);
          var index = result['redisIndex'];

          var fileUploadIndex = 0;

          // Remove the upload file from the redux store.
          for(var a = 0; a < fileUploads.length; ++a){

            if(fileUploads[a]['redisIndex'] == index){

              fileUploads.splice(a, 1);
              break;
            }
          }

          // Add the file upload result to the redux store.
          fileList.push(result);

          fileData['fileUploads'] = fileUploads;
          fileData['fileList'] = fileList;
          return fileData;

        case 'FILE_METADATA_RECEIVED':

          var fileData = Object.assign({}, state);

          for(var a = 0; a < action['fileList']['length']; ++a){

            var file = this.camelCaseIt(action['fileList'][a]);
            if(!action['fileList'][a].hasOwnProperty('finishTimestamp')){

              fileData['fileFails'].push(file);
            }
            else{

              fileData['fileList'].push(file);
            }
          }
          return fileData;

        case 'FILE_UPLOAD_PAUSED':

          var fileData = Object.assign({}, state);
          var fileUploads = fileData['fileUploads'];

          // MOD START
          var reqData = this.camelCaseIt(action['pauseResponse']['request']);
          for(var a = 0; a < fileUploads['length']; ++a){

            if(fileUploads[a]['redisIndex'] == reqData['redisIndex']){

              for(var key in reqData){

                fileUploads[a][key] = reqData[key];
              }
              break;
            }
          }
          return fileData;

        case 'QUEUE_UPLOAD':

          var fileData = Object.assign({}, state);
          var fileUploads = fileData['fileUploads'];

          for(var a = 0; a < fileUploads['length']; ++a){

            // Remove any 'queued' status
            if(fileUploads[a]['status'] == 'queued'){

              if(fileUploads[a]['currentBytePosition'] == 0){

                fileUploads[a]['status'] = 'initialized';
              }
              else{

                fileUploads[a]['status'] = 'paused';
              }
            }

            // Apply a 'queued' status to a matching file upload.
            if(fileUploads[a]['redisIndex'] == action['redisIndex']){

              fileUploads[a]['status'] = 'queued';
            }
          }
          return fileData;

        case 'FILE_UPLOAD_CANCELLED':

          var fileData = Object.assign({}, state);
          var fileUploads = fileData['fileUploads'];
          var cancelledFile = action['cancelledResponse']['request'];

          for(var a = 0; a < fileUploads['length']; ++a){

            if(fileUploads[a]['redisIndex'] == cancelledFile['redis_index']){

              fileUploads[a]['status'] = 'cancelled';
              fileData['fileFails'].push(fileUploads[a]);
              fileUploads.splice(a, 1);
            }
          }
          return fileData;

        case 'FILE_REMOVED':

          var fileData = Object.assign({}, state);
          var oldMetadata = action['oldMetadata'];

          for(var key in fileData){

            var fileRemoved = false;
            for(var a = 0; a < fileData[key]['length']; ++a){

              if(fileData[key][a]['redisIndex'] == oldMetadata['redis_index']){

                fileData[key].splice(a, 1);
                fileRemoved = true;
                break;
              }
            }

            if(fileRemoved){

              break;
            }
          }
          return fileData;

        default:

          var fileData = Object.assign({}, state);
          return fileData;
      }
    };
  }

  camelCaseIt(object){

    for(var key in object){

      object[CAMEL_CASE_IT(key)] = object[key];
      if(key.indexOf('_') != -1) delete object[key];
    }

    object = PARSE_REQUEST(object);
    return object;
  }
}