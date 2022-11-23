export type Attribute = {
  attribute_name: string;
  attribute_type: string;
  name: string;
  display_name: string;
  hidden: boolean;
  read_only: boolean;
  validation: Validation | null;
  restricted: boolean;
  model_name?: string;
  link_model_name?: string;
  description?: string;
};

export interface Validation {
  type: string;
  value: string | string[];
}

export interface Template {
  name: string;
  identifier: string;
  parent: string;
  attributes: {[key: string]: Attribute};
}

export interface Model {
  documents: {[key: string]: any};
  revisions: {[key: string]: any};
  template: Template;
  views: {[key: string]: any};
}
