import {
  checkStatus,
  handleFetchError,
  handleFetchSuccess,
  headers
} from 'etna-js/utils/fetch';

export type Attribute = {
  name: string;
  description?: string;
  attribute_name: string;
  attribute_type: string;
  attribute_group?: string;
  display_name?: string;
  restricted?: boolean;
  read_only?: boolean;
  hidden?: boolean;
  validation?: {[key: string]: any} | null;
};

export type LinkAttribute = Attribute & {
  link_model_name: string;
  link_attribute_name: string;
};

type AddAttributeParams = {
  model_name: string;
  attribute_name: string;
  description?: string;
  type: string;
  attribute_group?: string;
};

type AddAttributeAction = AddAttributeParams & {
  action_name: 'add_attribute';
};

type UpdateAttributeParams = {
  model_name: string;
  attribute_name: string;
  new_attribute_name?: string;
  description?: string;
  attribute_group?: string;
  display_name?: string;
  validation?: string;
  hidden?: boolean;
  read_only?: boolean;
  restricted?: boolean;
};

type UpdateAttributeAction = UpdateAttributeParams & {
  action_name: 'update_attribute';
};

type RenameAttributeAction = {
  action_name: 'rename_attribute';
  model_name: string;
  attribute_name: string;
  new_attribute_name: string;
};

type RemoveAttributeParams = {
  model_name: string;
  attribute_name: string;
};

type RemoveAttributeAction = RemoveAttributeParams & {
  action_name: 'remove_attribute';
};

// Basically params returned by the server but not
//   accepted as part of an update_attribute action.
const uneditableAttributes = [
  'desc',
  'name',
  'attribute_type',
  'options',
  'match'
];

const cleanPayload = (
  payload: (
    | AddAttributeAction
    | UpdateAttributeAction
    | RenameAttributeAction
    | RemoveAttributeAction
  )[]
) => {
  return payload.map((action: any) => {
    let clone = {...action};

    uneditableAttributes.forEach((attr) => {
      delete clone[attr];
    });

    return clone;
  });
};

const updateModel = (
  payload: (
    | AddAttributeAction
    | UpdateAttributeAction
    | RenameAttributeAction
    | RemoveAttributeAction
  )[]
) => {
  payload = cleanPayload(payload);
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

export const addAttribute = (params: AddAttributeParams) => {
  return updateModel([
    {
      action_name: 'add_attribute',
      ...params
    }
  ]);
};

export const updateAttribute = (params: UpdateAttributeParams) => {
  let updateAction: UpdateAttributeAction = {
    action_name: 'update_attribute',
    ...params
  };
  delete updateAction.new_attribute_name;

  let actions: (UpdateAttributeAction | RenameAttributeAction)[] = [
    updateAction
  ];

  if (
    params.new_attribute_name &&
    params.new_attribute_name != params.attribute_name
  ) {
    actions.push({
      action_name: 'rename_attribute',
      model_name: params.model_name,
      attribute_name: params.attribute_name,
      new_attribute_name: params.new_attribute_name
    });
  }

  return updateModel(actions);
};

export const removeAttribute = (params: RemoveAttributeParams) => {
  let removeAction: RemoveAttributeAction = {
    action_name: 'remove_attribute',
    ...params
  };

  let actions: RemoveAttributeAction[] = [removeAction];

  return updateModel(actions);
};
