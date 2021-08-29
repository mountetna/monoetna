declare module 'json-source-map';

type Job = {
  name: string;
  schema: any
};

type Etl = {
  name: string,
  etl: string,
  config: any,
  ran_at: string,
  run: number,
  output: string,
  status: string
};

