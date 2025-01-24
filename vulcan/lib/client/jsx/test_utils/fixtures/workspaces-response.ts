import {VulcanConfig, VulcanConfigElement, VulcanConfigRaw, WorkspacesResponse, WorkspaceStep} from '../../api_types';

export const test_step: WorkspaceStep = {
    name: 'count',
    input: {files: ['poem1', 'poem2']},
    output: {files: ['poem_count', 'poem_count_2']},
};

export const test_step_ui: WorkspaceStep = {
    name: 'count_selector',
    input: {files: ['poem1', 'poem2']},
    output: {files: ['poem_count', 'poem_count_2']},
    label: 'Count Type Selection',
    ui_component: 'dropdown',
    doc: 'help doc',
}

export const vc_test_step_ui: VulcanConfigElement = {
    name: 'count_selector',
    input: {files: ['options.txt']},
    output: {files: ['method']},
    display: 'Count Type Selection',
    ui_component: 'dropdown',
    doc: 'help doc',
}

export const vc_test_step_output: VulcanConfigElement = {
    name: 'count_display',
    output: {files: ['poem_count', 'poem_count_2']},
    display: 'Results',
    ui_component: 'raw',
}

export const test_vulcan_config_raw: VulcanConfigRaw = [
    vc_test_step_ui,
    vc_test_step_output,
]

export const test_vulcan_config: VulcanConfig = {
    count_selector: vc_test_step_ui,
    count_display: vc_test_step_output,
}

export const workspacesResponse: WorkspacesResponse = {
    workspaces: [
        {
            workspace_id: 1,
            workflow_id: 1,
            project: 'labors',
            steps: {
                select_prep: {
                    name: 'select_prep',
                    input: {},
                    output: {files: ['options.txt']},
                },
                count_selector: test_step_ui,
                count: test_step
            },
            dag: ['select_prep', 'count_selector', 'count'],
            last_config: {},
            last_job_status: {},
            author: 'someone',
            tags: ['public'],
            vulcan_config: test_vulcan_config
        }
    ]
};

