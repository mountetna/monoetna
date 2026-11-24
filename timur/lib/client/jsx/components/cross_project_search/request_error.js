import {json_error} from 'etna-js/utils/fetch';

export const REQUEST_FAILED = 'Request failed';

const oneLine = (message) => {
  if (typeof message === 'string') return message || REQUEST_FAILED;
  if (Array.isArray(message)) return message.join(', ') || REQUEST_FAILED;
  return REQUEST_FAILED;
};

// etna-js fetches reject with a promise of the response body; anything else
// (a network failure, a thrown Error) is not a promise. Resolves to one line of text.
export const requestError = (error) => {
  if (!error || typeof error.then !== 'function') {
    return Promise.resolve((error && error.message) || REQUEST_FAILED);
  }
  return new Promise((resolve) => {
    json_error((message) => resolve(oneLine(message)))(error)
      .catch(() => resolve(REQUEST_FAILED));
  });
};
