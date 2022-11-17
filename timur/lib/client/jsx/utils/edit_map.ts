import {Model, Attribute} from '../api/magma_api';

export const SNAKE_CASE = '^[a-z][a-z0-9]*(_[a-z0-9]*)*$';

export const SNAKE_CASE_STRICT = '^[a-z]*(_[a-z]*)*$';

export const COMMA_SEP = '^[a-zA-Z0-9]*(,[a-zA-Z0-9]*)*$';

export const VALIDATION_TYPES = ['Array', 'Regexp'];

export const EDITABLE_ATTRIBUTE_TYPES = [
  'string',
  'date_time',
  'shifted_date_time',
  'boolean',
  'float',
  'integer',
  'file',
  'image',
  'file_collection'
];

const CHILD_ATTRIBUTE_TYPES = ['table', 'child', 'collection'];

export const isLeafModel = (model: Model) => {
  return (
    Object.values(model.template.attributes)
      .map((attribute: Attribute): string => attribute.attribute_type)
      .filter((attr_type) => CHILD_ATTRIBUTE_TYPES.includes(attr_type))
      .length === 0
  );
};
