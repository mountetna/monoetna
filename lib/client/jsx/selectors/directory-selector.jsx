import { fileKey } from '../utils/file';

export const selectFiles = ({ directory: {files} }) => files;
export const selectBuckets = ({ directory: {buckets} }) => buckets;
export const selectFolders = ({ directory: {folders} }) => folders;
export const selectUploads = ({ directory: {uploads} }) => uploads;

export const selectUpload = ({ directory: {uploads} }, upload) => uploads[fileKey(upload)]
