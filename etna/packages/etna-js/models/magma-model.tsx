export type Model = {
  documents: {};
  template: {attributes: {[key: string]: Attribute}};
};

export type Attribute = {
  attribute_name: string;
  attribute_type: string;
  name: string;
  display_name: string;
  hidden: boolean;
  read_only: boolean;
  validation: {
    value: any;
  } | null;
  restricted: boolean;
  model_name?: string;
  link_model_name?: string;
};
