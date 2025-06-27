import {json_get, json_delete} from 'etna-js/utils/fetch';

export const addUsersSensor = (setInfo: Function) => json_get(
  `${CONFIG.janus_host}/api/admin/${CONFIG.project_name}/info`
).then(
  ({project}) => {
    const summary = project.permissions.map( ({role}:{role: string}) => role ).reduce(
      (s: any,r: string) => {
        s[r] = (s[r] || 0) + 1;
        s.count = (s.count || 0) + 1;
        return s
      }, {}
    );
    const count = (c: number,r: string) => c == 1 ? `1 ${r}` : `${c} ${r}s`;
    const level = (summary.count < 3) ? (summary.count < 2 ? 0 : 1) : 2;
    const text = [ 'administrator', 'editor', 'viewer', 'guest' ].map(
      role => role in summary ? count(summary[role], role) : ''
    ).filter(_=>_).join(', ')
    setInfo({level, text})
  }
) 


export const editModelsSensor = (setInfo: Function) => json_get(
  `${CONFIG.janus_host}/api/admin/${CONFIG.project_name}/info`
).then(
  ({project}) => {
    const summary = project.permissions.map( ({role}:{role: string}) => role ).reduce(
      (s: any,r: string) => {
        s[r] = (s[r] || 0) + 1;
        s.count = (s.count || 0) + 1;
        return s
      }, {}
    );
    const count = (c: number,r: string) => c == 1 ? `1 ${r}` : `${c} ${r}s`;
    const level = (summary.count < 3) ? (summary.count < 2 ? 0 : 1) : 2;
    const text = [ 'administrator', 'editor', 'viewer', 'guest' ].map(
      role => role in summary ? count(summary[role], role) : ''
    ).filter(_=>_).join(', ')
    setInfo({level, text})
  }
) 

export const editRulesSensor = (setInfo: Function) => setInfo(
  { level: 0, text: '1 of 12 models have identifier rules' }
);

export const createBucketsSensor = (setInfo: Function) => setInfo(
  { level: 2, text: '2 buckets created' }
);

export const addFilesSensor = (setInfo: Function) => setInfo(
  { level: 2, text: '572 files, 2.03 TB stored'}
);

export const linkRecordsSensor = (setInfo: Function) => setInfo(
      { level: 2, text: '200 records created' }
    );

export const createLoadersSensor = (setInfo: Function) => setInfo(
  { level: 2, text: '1 data loader, last run 2025-02-02' }
);

export const addWorkflowsSensor = (setInfo: Function) => setInfo(
  { level: 2, text: '1 workflow' }
);

export const runWorkflowsSensor = (setInfo: Function) => setInfo(
      { level: 2, text: '1 workspace, last run 2025-03-03' }
    );
