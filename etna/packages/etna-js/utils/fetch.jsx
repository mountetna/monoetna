export const checkStatus = (response) => {
  let content = isJSON(response) ? response.json() : response.text();
  if (response.status >= 200 && response.status < 300) {
    return content;
  } else {
    throw content;
  }
};

const isJSON = (response) =>
  response.headers.get('Content-Type') == 'application/json';

export const headers = (...types) => {
  var _headers = {};

  var add = (header, value) => (_headers[header] = value);

  for (var type of types) {
    switch (type) {
      case 'json':
        add('Content-Type', 'application/json');
        add('Accept', 'application/json');
        break;
      case 'csrf':
        let csrf = document.querySelector('meta[name=csrf-token]');
        if (csrf) add('X-CSRF-Token', csrf.getAttribute('content'));
        break;
    }
  }

  return _headers;
};

export const json_fetch = (method) => (path, params) =>
  fetch(path, {
    method,
    credentials: 'include',
    headers: headers('json'),
    ...(params && { body: JSON.stringify(params) })
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
