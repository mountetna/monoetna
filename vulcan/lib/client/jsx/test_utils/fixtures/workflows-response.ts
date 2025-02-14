import {WorkflowsResponse} from '../../api_types';

export const workflowsResponse: WorkflowsResponse = [
    {
        'id': 1,
        'name':'test_workflow',
        'projects': ['labors'],
        'tags':['demo', 'test'],
        'authors':['someone'],
        'branch': 'aaaa',
        'repo_remote_url': 'monoetna/test_workflow',
        'created_at': 1,
        'updated_at': 1,
        'displayName':'Test workflow'
    },
    {
        'id': 2,
        'name':'test_concurrent_workflow',
        'projects': ['labors'],
        'tags':['demo', 'test'],
        'authors':['someone'],
        'branch': 'aaaa',
        'repo_remote_url': 'monoetna/test_concurrent_workflow',
        'created_at': 2,
        'updated_at': 2,
        'displayName':'Test workflow 2'
    },
];