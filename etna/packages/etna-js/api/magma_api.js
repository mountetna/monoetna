import {checkStatus, headers} from '../utils/fetch';

export const magmaPath = (endpoint) => `${CONFIG.magma_host}/${endpoint}`;

const magmaPost = (endpoint, fetch, params) => {
  return fetch(magmaPath(endpoint), {
      method: 'POST',
      credentials: 'include',
      headers: headers('json'),
      body: JSON.stringify({
        project_name: CONFIG.project_name,
        ...params
      })
    })
    .then(checkStatus);
};

const create = (name, attributes) => {
  let element = document.createElement(name);
  for (let key in attributes) {
    element.setAttribute(key, attributes[key]);
  }
  return element;
};

const input = (name, value) => {
  return create('input', {type: 'hidden', name, value});
};

export const getTSVForm = ({
  model_name,
  filter,
  show_disconnected,
  attribute_names,
  output_predicate,
  expand_matrices,
  transpose
}) => {
  let {Authorization} = headers('auth');
  let data = {
    'X-Etna-Authorization': Authorization,
    project_name: CONFIG.project_name,
    model_name,
    record_names: 'all',
    attribute_names,
    filter,
    show_disconnected,
    format: 'tsv',
    output_predicate,
    expand_matrices,
    transpose
  };

  let form = create('form', {
    action: magmaPath('retrieve'),
    method: 'POST'
  });

  for (let name in data) {
    let value = data[name];
    if (value != undefined && value != null) {
      if (value instanceof Array) {
        value.forEach((val) => {
          form.appendChild(input(name + '[]', val));
        });
      } else {
        form.appendChild(input(name, value));
      }
    }
  }

  form.style.display = 'none';
  document.body.appendChild(form);
  form.submit();
  document.body.removeChild(form);
};

export const getDocuments = (doc_args, fetch=window.fetch) => {
  return magmaPost('retrieve', fetch, doc_args);
};

export const getModels = (project_name, fetch=window.fetch) => {
    return getDocuments(
      {
        project_name,
        model_name: 'all',
        record_names: [],
        attribute_names: 'all'
      },
      fetch
    );
}

export const postRevisions = (revision_data, fetch=window.fetch) => {
  return magmaPost('update', fetch, revision_data);
};

export const getAnswer = (question, fetch) => {
  return magmaPost('query', fetch, question);
};
