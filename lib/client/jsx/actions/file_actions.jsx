export const removeFile = () => (dispatch) => {
          if(action.fileMetadata == undefined) return;
          if(!confirm('Are you sure you want to remove this file?')) return;
          this.removeServerFiles('/remove-file', action.fileMetadata);
}

export const removeFailed = () => (dispatch) => {
   this.removeServerFiles('/remove-failed', action.fileMetadata);
}
