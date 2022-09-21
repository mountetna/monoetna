import { fileKey } from 'etna-js/utils/file';

const removeFiles = (old_files, action) => {
  let { files } = action;

  let removed_keys = files.reduce(
    (rk,file) => {
      rk[fileKey(file)] = true;
      return rk;
    }, {}
  );

  return Object.keys(old_files).filter(key => !(key in removed_keys)).reduce(
    (nf,key) => { nf[key] = old_files[key]; return nf; }, {}
  );
};

const addBuckets = (old_buckets, action) => {
  let { buckets } = action;

  let new_buckets = buckets.reduce((c, bucket) => {
    c[bucket.bucket_name] = bucket;
    return c;
  }, {});

  return {
    ...old_buckets,
    ...new_buckets
  };
};

const removeBucket = (old_buckets, action) => {
  let { bucket } = action;

  let { [bucket.bucket_name]: deleted_bucket, ...new_buckets } = old_buckets;

  return new_buckets;
};

const buckets = (old_buckets, action) => {
  if (!old_buckets) old_buckets = {};

  switch(action.type) {
    case 'REMOVE_BUCKET':
      return removeBucket(old_buckets, action);
    case 'ADD_BUCKETS':
      return addBuckets(old_buckets, action);
    default:
      return old_buckets;
  }
};

export default buckets;
