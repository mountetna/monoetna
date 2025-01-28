export type EtnaError = {
  error: string;
};

export type Job = {
  name: string;
  schema: any;
  runtime_params: any;
  secrets: any;
};

export type Workflow = {
  project_name: string;
  workflow_name: string;
  config_id: number;
  version_number: number;
  workflow_type: string;
  config: any;
  secrets: any;
  created_at: string;
};

export type Runtime = {
  config_id: number;
  config: any;
  run_interval: number;
};

export type Run = {
  run_id: number;
  created_at: string;
  config_id: number;
  version_number: number;
  status: string;
  finished_at: string;
};

export type Status = {
  workflow_type: string;
  workflow_name: string;
  config_id: number;
  pipeline_state: string;
  pipeline_started_at: string;
  pipeline_finished_at: string;
};

type WorkflowRevision = {
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
