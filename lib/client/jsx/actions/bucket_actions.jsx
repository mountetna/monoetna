import { postRetrieveBuckets, postUpdateBucket, postCreateBucket } from '../api/buckets_api';

const addBuckets = (buckets) => ({ type: 'ADD_BUCKETS', buckets });

export const retrieveBuckets = () => (dispatch) =>
  postRetrieveBuckets(CONFIG.project_name).then(
    ({buckets}) => dispatch(addBuckets(buckets))
  );

export const updateBucket = ({bucket}) => (dispatch) =>
  postUpdateBucket(CONFIG.project_name, bucket).then(
    ({bucket}) => dispatch(addBuckets([bucket]))
  );

export const createBucket = ({bucket}) => (dispatch) =>
  postCreateBucket(CONFIG.project_name, bucket).then(
    ({bucket}) => dispatch(addBuckets([bucket]))
  );
