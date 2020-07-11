import { fileKey } from '../utils/file';

export const selectFiles = ({ directory: {files} }, folder_name = null) => {
  if (folder_name == null) return files;
  const result = {};
  Object.keys(files).forEach(k => {
    let file = files[k];
    if (file.file_path === folder_name + '/' + file.file_name) {
      result[k] = file;
    }
  })
  return result;
};
export const selectBuckets = ({ directory: {buckets} }) => buckets;
export const selectFolders = ({ directory: {folders} }) => folders;
export const selectUploads = ({ directory: {uploads}
}) => uploads;

export const selectUpload = ({ directory: {uploads} }, upload) => uploads[fileKey(upload)]
