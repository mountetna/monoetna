import { postRetrieveBuckets, postUpdateBucket, postCreateBucket, deleteBucket } from '../api/buckets_api';
import { errorMessage } from './message_actions';

const addBuckets = (buckets) => ({ type: 'ADD_BUCKETS', buckets });
const removeBucket = (bucket) => ({ type: 'REMOVE_BUCKET', bucket });

export const retrieveBuckets = () => (dispatch) =>
  postRetrieveBuckets(CONFIG.project_name).then(
    ({buckets}) => dispatch(addBuckets(buckets))
  )
  .catch(
    errorMessage(dispatch, 'warning', 'Bucket updating failed', error => error)
  );


export const updateBucket = ({bucket}) => (dispatch) =>
  postUpdateBucket(CONFIG.project_name, bucket).then(
    ({bucket}) => dispatch(addBuckets([bucket]))
  )
  .catch(
    errorMessage(dispatch, 'warning', 'Bucket updating failed', error => error)
  );

export const createBucket = ({bucket}) => (dispatch) =>
  postCreateBucket(CONFIG.project_name, bucket).then(
    ({bucket}) => dispatch(addBuckets([bucket]))
  )
  .catch(
    errorMessage(dispatch, 'warning', 'Bucket creation failed', error => error)
  );

export const destroyBucket = ({bucket}) => (dispatch) => {
  if (!confirm(`Are you sure you want to remove ${bucket.bucket_name}?`)) return;

  deleteBucket(
    CONFIG.project_name, bucket.bucket_name
  )
    .then(({bucket}) => dispatch(removeBucket(bucket)))
    .catch(
      errorMessage(dispatch, 'warning', 'Bucket removal failed', error => error)
    );
};

