declare type ConfigType = {
  // projectName: string;
  token_name: string;
  janus_host: string;
  timur_host: string;
  metis_host: string;
  vulcan_host: string;
};

declare const CONFIG: ConfigType;

declare const ROUTES: {
  workflow(projectName: string, workflowName: string): string;
  fetch_workflows(): string;
  workflow_vignette(workflowName: string): string;
  submit(projectName: string, workflowName: string): string;
  status(projectName: string, workflowName: string): string;
  fetch_figures(projectName: string): string;
  create_figure(projectName: string): string;
  update_figure(projectName: string, figure_id: number): string;
  delete_figure(projectName: string, figure_id: number): string;
  fetch_figure(projectName: string, figure_id: number): string;
};
