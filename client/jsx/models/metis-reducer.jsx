export default class MetisReducer{

  reducer(){

    return (state = {}, action)=>{

      // Extract some widely use items.
      var fileData = Object.assign({}, state);
      var fileUploads = fileData['fileUploads'];
      var fileList = fileData['fileList'];
      var fileFails = fileData['fileFails'];

      switch(action['type']){

        case 'FILE_SELECTED':

          // Copy the selected file data to 'fileUploads' object.
          action['fileObject']['fileName'] = action['fileObject']['name'];
          action['fileObject']['originalName'] = action['fileObject']['name'];
          action['fileObject']['fileSize'] = action['fileObject']['size'];
          action['fileObject']['currentBytePosition'] = 0;
          action['fileObject']['status'] = 'unauthorized';
          action['fileObject']['reactKey'] = GENERATE_RAND_KEY();
          fileData['fileUploads'].push(action['fileObject']);
          break;

        case 'FILE_UPLOAD_AUTHORIZED':

          var response = this.camelCaseIt(action['authResponse']['request']);
          var index = this.getMatchingUploadIndex(fileUploads, response);

          // Append the HMAC signature and set the server current byte to 0.
          fileUploads[index]['hmacSignature'] = response['hmacSignature'];
          fileUploads[index]['currentBytePosition'] = 0;

          // Append all of the request items to the local file object.
          fileUploads[index] = Object.assign(fileUploads[index], response);
          break;

        case 'FILE_INITIALIZED':
        case 'FILE_UPLOAD_ACTIVE':
        case 'FILE_UPLOAD_PAUSED':

          var response = this.camelCaseIt(action['response']['request']);
          var index = this.getMatchingUploadIndex(fileUploads, response);

          // Append all of the request items to the local file object.
          fileUploads[index] = Object.assign(fileUploads[index], response);
          break;

        case 'FILE_UPLOAD_COMPLETE':

          var response =this.camelCaseIt(action['completeResponse']['request']);
          var index = this.getMatchingUploadIndex(fileUploads, response);

          /*
           * Move the completed upload metadata from the 'uploads' array to the 
           * list array.
           */
          fileUploads.splice(index, 1);
          fileList.push(response);
          break;

        case 'FILE_UPLOAD_CANCELLED':

          var response = this.camelCaseIt(action['cancelResponse']['request']);
          var index =this.getMatchingUploadIndex(fileUploads, response);

          /*
           * Move the cancelled upload metadata from the 'uploads' array to the 
           * failed array.
           */
          fileUploads.splice(index, 1);
          fileFails.push(cancelledResponse);
          break;

        case 'FILE_REMOVED':

          var response = this.camelCaseIt(action['response']);
          var index = this.getMatchingUploadIndex(fileUploads, response);

          // Remove the deleted item from the fileList
          fileList.splice(index, 1);
          break;

        case 'CLEAR_UPLOAD':

          for(var a = 0; a < fileUploads['length']; ++a){

            if(action['response']['reactKey'] == fileUploads[a]['reactKey']){

              fileUploads.splice(a, 1);
              break;
            }
          }
          break;

        case 'FILE_METADATA_RECEIVED':

          for(var a = 0; a < action['fileList']['length']; ++a){

            action['fileList'][a] = this.camelCaseIt(action['fileList'][a]);
            action['fileList'][a]['reactKey'] = GENERATE_RAND_KEY();

            if(!action['fileList'][a].hasOwnProperty('finishUpload')){

              fileData['fileFails'].push(action['fileList'][a]);
            }
            else{

              fileData['fileList'].push(action['fileList'][a]);
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
            if(fileUploads[a]['reactKey'] == action['reactKey']){

              fileUploads[a]['status'] = 'queued';
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