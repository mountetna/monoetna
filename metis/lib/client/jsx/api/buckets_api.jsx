import { json_get, json_delete, json_post } from 'etna-js/utils/fetch';

export const postRetrieveBuckets = (project_name) =>
  json_get(`/${project_name}/list`);

export const postUpdateBucket = (project_name, bucket) =>
  json_post(`/${project_name}/bucket/update/${bucket.bucket_name}`, bucket);

export const postCreateBucket = (project_name, bucket) =>
  json_post(`/${project_name}/bucket/create/${bucket.bucket_name}`, bucket);

export const deleteBucket = (project_name, bucket_name) =>
  json_delete(`/${project_name}/bucket/remove/${bucket_name}`);
