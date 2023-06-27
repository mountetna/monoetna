export type File = {
  file_name: string;
  file_path: string;
  bucket_name: string;
  project_name: string;
  read_only: boolean;
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
