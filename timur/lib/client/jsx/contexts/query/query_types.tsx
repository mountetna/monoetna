import {Attribute} from 'etna-js/models/magma-model';

export interface QuerySubclause {
  attributeName: string;
  attributeType: string;
  operator: string;
  operand: string | number;
}

export interface QueryClause {
  modelName: string;
  any: boolean;
  subclauses?: QuerySubclause[];
  // For migrating old, saved queries
  attributeName?: string;
  attributeType?: string;
  operator?: string;
  operand?: string | number;
}

export const EmptyQuerySubclause: QuerySubclause = {
  attributeName: '',
  attributeType: '',
  operator: '',
  operand: ''
};

export const EmptyQueryClause: QueryClause = {
  subclauses: [
    {
      ...EmptyQuerySubclause
    }
  ],
  modelName: '',
  any: true
};

export interface QueryBase {
  modelName: string;
}

export interface QuerySlice extends QueryBase {
  clause: QueryClause;
}

export interface QueryFilterAnyMap {
  [modelName: string]: boolean;
}

export interface QueryFilter extends QueryBase {
  anyMap: QueryFilterAnyMap;
  clauses: QueryClause[];
}

export interface QueryColumn {
  model_name: string;
  attribute_name: string;
  display_label: string;
  slices: QuerySlice[];
  predicate?: string;
}

export interface QueryResponse {
  answer: any[];
  format: any[];
  type: string;
}

export const EmptyQueryResponse: QueryResponse = {
  answer: [],
  format: [],
  type: 'none'
};

export interface QueryTableColumn {
  label: string;
  colId: string;
  modelName: string;
  attribute: Attribute;
  matrixHeadings: string[];
  predicate: string | undefined;
}

interface InputQueryMap {
  [key: string]: string;
}

export interface Workflow {
  queryAction?: string;
  name: string;
  displayName: string;
  icon?: string;
  inputQueryMap?: InputQueryMap;
}

export interface QueryPayload {
  query: string | any[];
  expand_matrices?: boolean;
  user_columns?: string[];
  transpose?: boolean;
  format?: string;
}

export interface CreateFigurePayload {
  title: string;
  workflow_name: string;
  inputs: {[key: string]: any};
}

export interface SavedQuery {
  id: number;
  created_at: string;
  comment: string;
  query: string;
  unpackedQuery?: string;
}
