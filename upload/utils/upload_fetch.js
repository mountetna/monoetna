// Some of this is duplicated code from etna-js/utils/fetch
//   But we need to sandbox the web-worker code from
//   the client-side code.

export const isJSON = (response) =>
  response.headers.get('Content-Type') == 'application/json';

export const checkStatus = (response) => {
  let content = isJSON(response) ? response.json() : response.text();
  if (response.status >= 200 && response.status < 300) {
    return content;
  } else {
    throw content;
  }
};

export const json_fetch = (method) => (path, params) =>
  fetch(path, {
    method,
    credentials: 'include',
    headers: {
      'Content-type': 'application/json'
    },
    ...(params && {body: JSON.stringify(params)})
  }).then(checkStatus);

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
