import downloadjs from 'downloadjs';
import * as Cookies from './cookies';

export const checkStatus = (response) => {
  let content = isJSON(response) ? response.json() : response.text();
  if (response.status >= 200 && response.status < 300) {
    return content;
  } else {
    throw content;
  }
};

export const parseJSON = (response) => {
  try {
    return response.json();
  } catch (err) {
    console.log(err);
    throw err;
  }
};

export const makeBlob = (response) => {
  return response.blob();
};

export const generateDownload = (filename) => {
  return (blob) => {
    return downloadjs(blob, filename, blob.type);
  };
};

const isJSON = (response) =>
  response.headers.get('Content-Type') == 'application/json';

export const headers = (...types) => {
  let _headers = {};
  let add = (header, value) => (_headers[header] = value);

  for (let type of types) {
    switch (type) {
      case 'json':
        add('Content-Type', 'application/json');
        break;
      case 'csrf':
        let csrf = document.querySelector('meta[name=csrf-token]');
        if (csrf) add('X-CSRF-Token', csrf.getAttribute('content'));
        break;
      case 'auth':
        let token = Cookies.getItem(CONFIG.token_name);
        add('Authorization', `Etna ${token}`);
        break;
      default:
        break;
    }
  }

  return _headers;
};


const opts = {
  credentials: 'same-origin',
  headers: headers('json', 'csrf')
};

export const getOpts = { method: 'GET', ...opts };
export const postOpts = (body) => ({ method: 'POST', body: JSON.stringify(body), ...opts });
export const deleteOpts = { method: 'DELETE', ...opts };


export const json_fetch = (method) => (path, params) =>
  fetch((CONFIG.baseURL || '') + path, {
    method,
    credentials: 'include',
    headers: headers('json'),
    ...(params && {body: JSON.stringify(params)})
  }).then(checkStatus);

export const json_get = json_fetch('GET');
export const json_delete = json_fetch('DELETE');
export const json_post = json_fetch('POST');

export const form_post = (path, params, performCheckStatus = true) => {
  let form = new FormData();

  Object.keys(params).forEach((key) => form.append(key, params[key]));

  return fetch(path, {
    method: 'POST',
    credentials: 'same-origin',
    body: form
  }).then(performCheckStatus ? checkStatus : (v) => v);
};

export const handleFetchError = (e) => {
  console.log(e);
  return e.then((response) => {
    if (!response) {
      return Promise.reject([`Something is amiss. ${e}`]);
    }

    let errStr = response.error
      ? response.error
      : response.errors.map((error) => `* ${error}`);
    errStr = [`### Our request was refused.\n\n${errStr}`];
    return Promise.reject(errStr);
  });
};

export const handleFetchSuccess = (response) => {
  if (response && typeof response === 'object' && 'error' in response) {
    console.log(response.error);
    return Promise.reject([`There was a ${response.type} error.`]);
  }
  return Promise.resolve(response);
};
