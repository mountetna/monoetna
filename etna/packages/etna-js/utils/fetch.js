import downloadjs from 'downloadjs';
import * as Cookies from './cookies';
import {base_json_fetch, checkStatus, isJSON, form_post} from './base_fetch';

export {checkStatus, isJSON, form_post};

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

export const getOpts = {method: 'GET', ...opts};
export const postOpts = (body) => ({
  method: 'POST',
  body: JSON.stringify(body),
  ...opts
});
export const deleteOpts = {method: 'DELETE', ...opts};

export const json_fetch = (method) => (path, params) =>
  base_json_fetch(method)(`${CONFIG.baseURL || ''}${path}`, params);

export const json_get = json_fetch('GET');
export const json_delete = json_fetch('DELETE');
export const json_post = json_fetch('POST');

export const json_error = handler => e => e.then(response => handler(error_string(response)));

const error_string = (response,md=false) => !response ? "Request error" : response.error
  ? response.error
  : response.errors
  ? response.errors.map(e => (md ? '* ' : '') + `${e.message ? e.message : e}`)
  : response;

export const handleFetchError = (e) => {
  console.error(e);
  return Promise.resolve(e).then((response) => {
    if (!response) {
      return Promise.reject([`Something is amiss. ${e}`]);
    }

    let errStr = error_string(response,true);
    errStr = [`### Our request was refused.\n\n${errStr}`];
    return Promise.reject(errStr);
  });
};

export const handleFetchSuccess = (response) => {
  if (response && typeof response === 'object' && 'error' in response) {
    console.error(response.error);
    return Promise.reject([
      `There was a ${response.type} error. ${response.error}`
    ]);
  }
  return Promise.resolve(response);
};
