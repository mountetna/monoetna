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

export type MetisFile = {
  id: number;
  file_name: string;
  file_path: string;
  bucket_name: string;
  project_name: string;
  read_only: boolean;
  updated_at: string;
  download_url: string;
};


export type Script = any;

export type ScriptItem = {
  type: string;
  projectName: string;
  bucketName: string;
  modelName: string;
  value: any;
  script: Script;
  update: Function;
  classes: any;
};
