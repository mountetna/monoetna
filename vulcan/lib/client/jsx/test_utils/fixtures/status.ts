import { WorkspaceStatus } from '../../api_types';
import { Maybe } from '../../selectors/maybe';

export const statusEmpty: WorkspaceStatus = {
    steps: {
        select_prep: {
            name: 'select_prep',
            status: 'pending',
            statusFine: 'NOT STARTED',
        },
        count_selector: {
            name: 'count_selector',
            status: 'pending',
            statusFine: 'NOT STARTED',
        },
        count: {
            name: 'count',
            status: 'pending',
            statusFine: 'NOT STARTED',
        },
    },
    output_files: [],
    file_contents: {},
    last_params: {},
    params: {},
    ui_contents: {},
};

export const statusCompleted: WorkspaceStatus = {
    steps: {
        select_prep: {
            name: 'select_prep',
            status: 'complete',
            statusFine: 'COMPLETED',
        },
        count_selector: {
            name: 'count_selector',
            status: 'complete',
            statusFine: 'COMPLETED',
        },
        count: {
            name: 'count',
            status: 'complete',
            statusFine: 'COMPLETED',
        },
    },
    output_files: ['options.txt', 'method', 'poem1', 'poem2', 'poem_count', 'poem_count_2'],
    file_contents: {
        'options.txt': '"bytes", "chars", "words"',
        'method': '"words"',
        poem_count: '4',
        poem_count_2: '6',
    },
    last_params: {},
    params: {},
    ui_contents: {
        count_selector: {
            method: ["words"] as Maybe<any>
        }
    },
}

export const statusPartiallyComplete: WorkspaceStatus = {
    steps: {
        select_prep: {
            name: 'select_prep',
            status: 'complete',
            statusFine: 'COMPLETED',
        },
        count_selector: {
            name: 'count_selector',
            status: 'pending',
            statusFine: 'NOT STARTED',
        },
        count: {
            name: 'count',
            status: 'pending',
            statusFine: 'NOT STARTED',
        },
    },
    output_files: ['options.txt', 'poem1', 'poem2'],
    file_contents: {
        'options.txt': '"bytes", "chars", "words"',
    },
    last_params: {},
    params: {},
    ui_contents: {
        count_selector: {
            method: null as Maybe<any>
        }
    },
}

export const statusWithError: WorkspaceStatus = {
    ...statusEmpty,
    steps: {
        ...statusEmpty['steps'],
        count: {
            name: 'count',
            status: 'pending',
            statusFine: 'NOT STARTED',
            error: "Python error while executing script: Traceback (most recent call last):\n  File \"/archimedes-exec/tmpwvdtot6a/script.py\", line 4, in <module>\n    b = int(open(input_path('b'), 'r').read())\nValueError: invalid literal for int() with base 10: 'abc'\n"
        }
    }
};
