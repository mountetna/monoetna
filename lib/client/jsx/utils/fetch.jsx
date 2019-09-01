export const checkStatus = (response) => {
  if (response.status >= 200 && response.status < 300) {
    return response
  } else {
    var error = new Error(response.statusText)
    error.response = response
    error.fetch = true;
    throw error
  }
}

export const parseJSON = (response) => response.json()

export const makeBlob = (response) => response.blob()

export const headers = (...types) => {
  var _headers = {}

  var add = (header, value) => _headers[header] = value

  for (var type of types) {
    switch(type) {
      case 'json':
        add( 'Content-Type', 'application/json');
        add( 'Accept', 'application/json');
        break;
      case 'csrf':
        let csrf = document.querySelector('meta[name=csrf-token]');
        if (csrf) add('X-CSRF-Token', csrf.getAttribute('content'));
        break;
    }
  }

  return _headers
}

export const json_fetch = (method) => (path, params) => fetch(path,
  {
    method,
    credentials: 'same-origin',
    headers: headers('json'),
    ...params && { body: JSON.stringify(params) }
  }).then(checkStatus).then(parseJSON);

export const json_get = json_fetch('GET');
export const json_delete = json_fetch('DELETE');
export const json_post = json_fetch('POST');

export const form_post = (path, params) => {
  let form = new FormData();

  Object.keys(params).forEach(key => form.append(key, params[key]));

  return fetch(path,
    {
      method: 'POST',
      credentials: 'same-origin',
      body: form
    }).then(checkStatus).then(parseJSON);
}
