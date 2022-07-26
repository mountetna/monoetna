declare module 'json-source-map';

type EtnaError = {
  error: string
};

type Job = {
  name: string,
  schema: any,
  params: any,
  secrets: any
};

type Etl = {
  project_name: string,
  name: string,
  etl: string,
  config: any,
  ran_at: string,
  updated_at: string,
  run_interval: number,
  output: string,
  secrets: any,
  params: any,
  archived: boolean,
  status: string
};

type EtlRevision = {
  config: any,
  comment: string,
  updated_at: string
}

type MagmaModel = any;

type MagmaModels = {
  [model_name: string]: MagmaModel
}
