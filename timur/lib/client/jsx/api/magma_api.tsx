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
  format_hint?: string;
  restricted?: boolean;
  read_only?: boolean;
  hidden?: boolean;
  validation?: {[key: string]: any} | null;
};

export type Model = {
  template: Template;
};

export type Template = {
  attributes: {
    [key: string]: Attribute
  };
}

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
  format_hint?: string;
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

type AddLinkParams = {
  modelName: string;
  linkAttributeName: string;
  reciprocalAttributeName: string;
  reciprocalModelName: string;
  reciprocalLinkType: 'link' | 'child' | 'collection';
};

type LinkAction = {
  model_name: string;
  attribute_name: string;
  type: 'link' | 'child' | 'collection';
};

type AddLinkAction = {
  action_name: 'add_link';
  links: [LinkAction, LinkAction];
};

type RemoveLinkParams = {
  model_name: string;
  attribute_name: string;
};

type RemoveLinkAction = RemoveLinkParams & {
  action_name: 'remove_link';
};

type AddModelParams = {
  model_name: string;
  identifier: string;
  parent_model_name: string;
  parent_link_type: 'child' | 'collection';
};

type AddModelAction = AddModelParams & {
  action_name: 'add_model';
};

type RemoveModelParams = {
  model_name: string;
};

type RemoveModelAction = RemoveModelParams & {
  action_name: 'remove_model';
};

type ReparentModelParams = {
  model_name: string;
  parent_model_name: string;
};

type ReparentModelAction = ReparentModelParams & {
  action_name: 'reparent_model';
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
    | AddLinkAction
    | RemoveLinkAction
    | AddModelAction
    | RemoveModelAction
    | ReparentModelAction
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
    | AddLinkAction
    | RemoveLinkAction
    | AddModelAction
    | RemoveModelAction
    | ReparentModelAction
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

export const addAttributes = (attributes: AddAttributeParams[]) => {
  return updateModel(attributes.map(params => (
    {
      action_name: 'add_attribute',
      ...params
    }
  )));
};

export const addAttribute = (params: AddAttributeParams) => addAttributes([params]);

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

export const addLink = (params: AddLinkParams) => {
  let addLinkAction: AddLinkAction = {
    action_name: 'add_link',
    links: [
      {
        model_name: params.modelName,
        attribute_name: params.linkAttributeName,
        type: 'link'
      },
      {
        model_name: params.reciprocalModelName,
        attribute_name: params.reciprocalAttributeName,
        type: params.reciprocalLinkType
      }
    ]
  };

  let actions: AddLinkAction[] = [addLinkAction];

  return updateModel(actions);
};

export const removeLink = (params: RemoveLinkParams) => {
  let removeAction: RemoveLinkAction = {
    action_name: 'remove_link',
    ...params
  };

  let actions: RemoveLinkAction[] = [removeAction];

  return updateModel(actions);
};

export const addModel = (params: AddModelParams) => {
  let addModelAction: AddModelAction = {
    action_name: 'add_model',
    ...params
  };

  let actions: AddModelAction[] = [addModelAction];

  return updateModel(actions);
};

export const removeModel = (params: RemoveModelParams) => {
  let removeModelAction: RemoveModelAction = {
    action_name: 'remove_model',
    ...params
  };

  let actions: RemoveModelAction[] = [removeModelAction];

  return updateModel(actions);
};

export const reparentModel = (params: ReparentModelParams) => {
  let reparentModelAction: ReparentModelAction = {
    action_name: 'reparent_model',
    ...params
  };

  let actions: ReparentModelAction[] = [reparentModelAction];

  return updateModel(actions);
};
