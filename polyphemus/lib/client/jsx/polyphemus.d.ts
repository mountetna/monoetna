export type EtnaError = {
  error: string;
};

export type Job = {
  name: string;
  schema: any;
  params: any;
  secrets: any;
};

export type Etl = {
  project_name: string;
  name: string;
  config_id: number;
  etl: string;
  config: any;
  ran_at: string;
  updated_at: string;
  run_interval: number;
  output: string;
  secrets: any;
  params: any;
  status: string;
};

type EtlRevision = {
  config: any;
  comment: string;
  updated_at: string;
};

type MagmaModel = any;

type MagmaModels = {
  [model_name: string]: MagmaModel;
};
