import {
  checkStatus,
  handleFetchError,
  handleFetchSuccess,
  headers
} from 'etna-js/utils/fetch';

type AddAttributeAction = {
  action_name: 'add_attribute';
  model_name: string;
  description?: string;
  attribute_name: string;
  type: string;
  attribute_group?: string;
};

type AddAttributeRequest = {
  model_name: string;
  name: string;
  description: string;
  type: string;
  group: string;
};

const updateModel = (payload: AddAttributeAction[]) => {
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

export const addAttribute = ({
  model_name,
  name,
  description,
  type,
  group
}: AddAttributeRequest) => {
  return updateModel([
    {
      action_name: 'add_attribute',
      model_name,
      description,
      attribute_name: name,
      type,
      attribute_group: group
    }
  ]);
};
