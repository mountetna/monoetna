import { headers, parseJSON, checkStatus,
  json_get, json_delete, json_post } from '../utils/fetch';

export const postRetrieveBuckets = (project_name) =>
  json_get(`/${project_name}/list`);

export const postUpdateBucket = (project_name, bucket) =>
  json_post(`/${project_name}/bucket/update/${bucket.bucket_name}`, bucket);

export const postCreateBucket = (project_name, bucket) =>
  json_post(`/${project_name}/bucket/create/${bucket.bucket_name}`, bucket);
