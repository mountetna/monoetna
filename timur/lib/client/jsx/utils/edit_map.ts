import {ModelObject, AttributeObject} from 'etna-js/models/magma-model';

export const SNAKE_CASE = '^[a-z][a-z0-9]*(_[a-z0-9]*)*$';

export const SNAKE_CASE_STRICT = '^[a-z]*(_[a-z]*)*$';

export const COMMA_SEP = '^[a-zA-Z0-9]*(,[a-zA-Z0-9]*)*$';

export const COMMA_SEP_WITH_SPACES = '^[a-zA-Z0-9 ]*(,[a-zA-Z0-9 ]*)*$';

export const VALIDATION_TYPES = ['Array', 'Regexp'];

export const EDITABLE_OPTIONS = [
    'description',
    'display_name',
    'format_hint',
    'hidden',
    'attribute_group',
    'index',
    'link_model_name',
    'read_only',
    'restricted',
    'unique',
    'validation',
    'link_attribute_name'
]

export const REMOVABLE_ATTRIBUTE_TYPES = [
  'string',
  'date_time',
  'shifted_date_time',
  'boolean',
  'float',
  'integer',
  'file',
  'image',
  'file_collection',
  'matrix'
];

export const EDITABLE_ATTRIBUTE_TYPES = [
  ...REMOVABLE_ATTRIBUTE_TYPES,
  'identifier'
];

export const UNREMOVABLE_ATTRIBUTE_NAMES = ['created_at', 'updated_at'];

export const UNEDITABLE_ATTRIBUTE_NAMES = ['created_at', 'updated_at'];

export const CHILD_ATTRIBUTE_TYPES = ['table', 'child', 'collection'];

export const isLeafModel = (model: ModelObject) => {
  return (
    Object.values(model.template.attributes)
      .map((attribute: AttributeObject): string => attribute.attribute_type)
      .filter((attr_type) => CHILD_ATTRIBUTE_TYPES.includes(attr_type))
      .length === 0
  );
};
