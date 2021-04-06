declare type ConfigType = {
    // project_name: string;
    token_name: string;
    janus_host: string;
    timur_host: string;
    metis_host: string;
    vulcan_host: string;
}

declare const CONFIG: ConfigType;

declare const ROUTES: {
    workflow(): string,
    fetch_workflows(): string,
    workflow_vignette(workflowName: string): string,
    submit(workflowName: string): string,
    status(workflowName: string): string,
}