export default class MetisReducer{

  reducer(){

    return (state = {}, action)=>{

      switch(action['type']){

        case 'FILE_SELECTED':
          
          var nextState = Object.assign({}, state);

          // MOD START
          var fileObject = action['data'];
          fileObject['file_name'] = fileObject['name'];
          fileObject['original_name'] = fileObject['name'];
          fileObject['file_size'] = fileObject['size'];
          fileObject['user_email'] = state['userInfo']['user_email'];
          fileObject['status'] = 'queued';

          /*
           * This item will be replaced by a Redis id. We will use it in Redis
           * as a unique entry identifier and also in React as a DOM key.
           * For new uploads we just generate a psudo random string, here at the
           * client, to hold a place.
           */
          fileObject['redis_index'] = GENERATE_RAND_KEY(); 
          // MOD END
          
          nextState['fileUploads'].push(fileObject);
          return nextState;

        case 'FILE_UPLOAD_AUTHORIZED':

          var fileUploads = state['fileUploads'];
          var nextState = Object.assign({}, state);

          // MOD START
          var authResponse = action['data'];
          var fileUpload = null;
          var fileUploadIndex = 0;

          // Select the file to upload from the redux store using the old index.
          var oldIndex = authResponse['request']['old_index'];
          for(var a = 0; a < fileUploads.length; ++a){

            if(fileUploads[a]['redis_index'] == oldIndex){

              fileUploadIndex = a; 
              fileUpload = fileUploads[a];
              break;
            }
          }
          delete authResponse['request']['old_index'];

          // Append the signature and set the server current byte to 0
          fileUpload['signature'] = authResponse['signature'];
          fileUpload['current_byte_position'] = 0;

          // Append all of the request items to the local file object
          for(var key in authResponse['request']){

            fileUpload[key] = authResponse['request'][key];
          }
          // MOD END

          nextState['fileUploads'][fileUploadIndex] = fileUpload;
          return nextState;

        case 'FILE_UPLOAD_ACTIVE':

          var fileUploads = state['fileUploads'];
          var nextState = Object.assign({}, state);

          // MOD START
          var response = action['data'];
          var index = response['request']['redis_index'];
          var fileUpload = null;
          var fileUploadIndex = 0;

          // Select the file to upload from the redux store.
          for(var a = 0; a < fileUploads.length; ++a){

            if(fileUploads[a]['redis_index'] == index){

              fileUploadIndex = a; 
              fileUpload = fileUploads[a];
              break;
            }
          }

          // Append the status to the file upload object.
          fileUpload['status'] = response['status'];

          // Append all of the request items to the local file object
          for(var key in response['request']){

            fileUpload[key] = response['request'][key];
          }

          nextState['fileUploads'][fileUploadIndex] = fileUpload;
          // MOD END

          return nextState;

        case 'FILE_UPLOAD_COMPLETE':

          var fileUploads = state['fileUploads'];
          var fileList = state['fileList'];
          var nextState = Object.assign({}, state);

          // MOD START
          var result = action['data']['result'];
          var index = result['redis_index'];

          var fileUploadIndex = 0;

          // Remove the upload file from the redux store.
          for(var a = 0; a < fileUploads.length; ++a){

            if(fileUploads[a]['redis_index'] == index){

              fileUploads.splice(a, 1);
              break;
            }
          }

          // Add the file upload result to the redux store.
          fileList.push(result);

          nextState['fileUploads'] = fileUploads;
          nextState['fileList'] = fileList;

          return nextState;
        default:

          var nextState = Object.assign({}, state);
          return nextState;
      }
    };
  }
}