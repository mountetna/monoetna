export default class MetisReducer{

  reducer(){

    return (state = {}, action)=>{

      var fileData = Object.assign({}, state);
      if('fileUploads' in fileData) var fileUploads = fileData['fileUploads'];

      switch(action['type']){

        case 'FILE_SELECTED':

          // Copy the selected file data to 'fileUploads' object.
          var fileObject = action['fileObject'];
          fileObject['fileName'] = fileObject['name'];
          fileObject['originalName'] = fileObject['name'];
          fileObject['fileSize'] = fileObject['size'];
          fileObject['currentBytePosition'] = 0;
          fileObject['status'] = 'unauthorized';
          fileObject['reactKey'] = GENERATE_RAND_KEY();
          fileData['fileUploads'].push(fileObject);
          break;

        case 'FILE_UPLOAD_AUTHORIZED':

          var authResponse = Object.assign({}, action['authResponse']);
          authResponse = this.camelCaseIt(authResponse['request']);

          // Find the local File Object.
          var index = this.getMatchingUploadIndex(fileUploads, authResponse);

          // Append the HMAC signature and set the server current byte to 0.
          fileUploads[index]['hmacSignature'] = authResponse['hmacSignature'];
          fileUploads[index]['currentBytePosition'] = 0;

          // Append all of the request items to the local file object.
          fileUploads[index] = Object.assign(fileUploads[index], authResponse);
          break;

        case 'FILE_INITIALIZED':

          var initResponse = Object.assign({}, action['initResponse']);
          initResponse = this.camelCaseIt(initResponse['request']);

          // Find the local File Object.
          var index = this.getMatchingUploadIndex(fileUploads, initResponse);

          // Append all of the request items to the local file object.
          fileUploads[index] = Object.assign(fileUploads[index], initResponse);
          break;

        case 'FILE_UPLOAD_ACTIVE':

          var fileUploads = fileData['fileUploads'];

          var response = action['uploadResponse'];
          var index = response['request']['redis_index'];
          var fileUpload = null;
          var fileUploadIndex = 0;

          // Select the file to upload from the redux store.
          for(var a = 0; a < fileUploads.length; ++a){

            if(fileUploads[a]['dbIndex'] == index){

              fileUploadIndex = a; 
              fileUpload = fileUploads[a];
              break;
            }
          }

          // Append all of the request items to the local file object.
          response['request'] = this.camelCaseIt(response['request']);
          fileUpload = Object.assign(fileUpload, response['request']);

          fileData['fileUploads'][fileUploadIndex] = fileUpload;
          break;

        case 'FILE_UPLOAD_COMPLETE':

          var fileList = fileData['fileList'];

          // MOD START
          var result = this.camelCaseIt(action['uploadResponse']['result']);
          var index = result['dbIndex'];

          var fileUploadIndex = 0;

          // Remove the upload file from the redux store.
          for(var a = 0; a < fileUploads.length; ++a){

            if(fileUploads[a]['dbIndex'] == index){

              fileUploads.splice(a, 1);
              break;
            }
          }

          // Add the file upload result to the redux store.
          fileList.push(result);

          fileData['fileUploads'] = fileUploads;
          fileData['fileList'] = fileList;
          break;

        case 'FILE_METADATA_RECEIVED':

          for(var a = 0; a < action['fileList']['length']; ++a){

            var file = this.camelCaseIt(action['fileList'][a]);
            file['reactKey'] = GENERATE_RAND_KEY();

            if(!action['fileList'][a].hasOwnProperty('finishTimestamp')){

              fileData['fileFails'].push(file);
            }
            else{

              fileData['fileList'].push(file);
            }
          }
          break;

        case 'FILE_UPLOAD_PAUSED':

          // MOD START
          var reqData = this.camelCaseIt(action['pauseResponse']['request']);
          for(var a = 0; a < fileUploads['length']; ++a){

            if(fileUploads[a]['dbIndex'] == reqData['dbIndex']){

              for(var key in reqData){

                fileUploads[a][key] = reqData[key];
              }
              break;
            }
          }
          break;

        case 'QUEUE_UPLOAD':

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
            if(fileUploads[a]['dbIndex'] == action['dbIndex']){

              fileUploads[a]['status'] = 'queued';
            }
          }
          break;

        case 'FILE_UPLOAD_CANCELLED':

          var cancelledFile = action['cancelledResponse']['request'];

          for(var a = 0; a < fileUploads['length']; ++a){

            if(fileUploads[a]['dbIndex'] == cancelledFile['redis_index']){

              fileUploads[a]['status'] = 'cancelled';
              fileData['fileFails'].push(fileUploads[a]);
              fileUploads.splice(a, 1);
            }
          }
          break;

        case 'FILE_REMOVED':

          var oldMetadata = action['oldMetadata'];

          for(var key in fileData){

            var fileRemoved = false;
            for(var a = 0; a < fileData[key]['length']; ++a){

              if(fileData[key][a]['dbIndex'] == oldMetadata['redis_index']){

                fileData[key].splice(a, 1);
                fileRemoved = true;
                break;
              }
            }

            if(fileRemoved){

              break;
            }
          }

          break;
        default:

          break;
      }

      return fileData;
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

  // Find the local File Object.
  getMatchingUploadIndex(fileUploads, responseData){

    var index = 0;
    for(var a = 0; a < fileUploads['length']; ++a){

      if((fileUploads[a]['fileName'] == responseData['fileName']) &&
        (fileUploads[a]['projectName'] == responseData['projectName']) &&
        (fileUploads[a]['groupName'] == responseData['groupName'])){

        index = a; 
        break;
      }
    }
    return index;
  }
}