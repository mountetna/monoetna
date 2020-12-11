import {fileKey} from '../utils/file';

export const selectFiles = ({directory: {files}}, folder_name = null) => {
  if (folder_name == null) return files;
  const result = {};
  Object.keys(files).forEach((k) => {
    let file = files[k];
    if (file.file_path === folder_name + '/' + file.file_name) {
      result[k] = file;
    }
  });
  return result;
};
export const selectBuckets = ({directory: {buckets}}) => buckets;
export const selectFolders = ({directory: {folders}}) => folders;
export const selectUploads = ({directory: {uploads}}) => uploads;
export const selectUpload = ({directory: {uploads}}, upload) =>
  uploads[fileKey(upload)];

export const selectUploadForRevision = (
  uploads,
  model_name,
  record_name,
  attribute_name
) => {
  let match = null;
  Object.keys(uploads).forEach((file_key) => {
    let upload = uploads[file_key];
    if (
      upload.model_name === model_name &&
      upload.record_name === record_name &&
      upload.attribute_name === attribute_name
    ) {
      match = upload;
    }
  });
  return match;
};
