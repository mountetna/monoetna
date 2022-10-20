import {
  checkStatus,
  handleFetchError,
  handleFetchSuccess,
  headers
} from 'etna-js/utils/fetch';

export const addAttribute = ({model_name, name, description, type}) => {
  return updateModel([
    {
      action_name: 'add_attribute',
      model_name,
      description,
      attribute_name: name,
      type
    }
  ]);
};

const updateModel = (payload: any) => {
  return fetch(`${CONFIG.magma_host}/update_model`, {
    method: 'POST',
    credentials: 'include',
    headers: headers('json'),
    body: JSON.stringify({
      actions: payload,
      project_name: CONFIG.project_name
    })
  })
    .then(checkStatus)
    .then(handleFetchSuccess)
    .catch(handleFetchError);
};
