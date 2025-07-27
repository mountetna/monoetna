import {json_get, json_delete} from 'etna-js/utils/fetch';
import {byteFormat} from 'etna-js/utils/format';

const plural = (word,count) => count == 1 ? word : word + 's';
export const addUsersSensor = (setInfo: Function, state: any) => {
  const { projectInfo } = state;

  if (!projectInfo) return;

  const summary = projectInfo.permissions.map( ({role}:{role: string}) => role ).reduce(
    (s: any,r: string) => {
      s[r] = (s[r] || 0) + 1;
      s.count = (s.count || 0) + 1;
      return s
    }, {}
  );
  const level = (summary.count < 3) ? (summary.count < 2 ? 0 : 1) : 2;
  const text = [ 'administrator', 'editor', 'viewer', 'guest' ].map(
    role => role in summary ? `${summary[role]} ${plural(role, summary[role])}` : ''
  ).filter(_=>_).join(', ');

  setInfo({level, text})
}

export const editModelsSensor = (setInfo: Function, state: any) => {
  const { models } = state;

  if (!models) return;
  console.log({models});

  const modelCount = Object.values(models).length;
  const filledCount = Object.values(models).filter(
    model => Object.values(model.template.attributes).filter(
      attribute => !attribute.link_model_name && attribute.attribute_type != 'identifier' && attribute.attribute_name != 'created_at' && attribute.attribute_name != 'updated_at'
    ).length > 0
  ).length

  setInfo({
    level: modelCount < 3 ? 0 : Math.min(Math.floor((filledCount / modelCount) * 3), 2),
    text: `${filledCount} of ${modelCount} ${plural('model',modelCount)} with attributes`
  });
} 

export const editRulesSensor = (setInfo: Function, state: any) => {
  const { models, rules } = state;

  if (!models || !rules) return;
  console.log({rules});

  const modelsWithIds = Object.keys(models).filter(
    modelName => {
      const template = models[modelName].template;
      const id = template.attributes[template.identifier];
      return (id.attribute_name != 'id' && !id.hidden)
    }
  );

  const modelsWithNames = modelsWithIds.filter(
    modelName => modelName in rules
  );
    
  setInfo({
    level: modelsWithIds.length == 0 ? 0 : Math.min(Math.floor(modelsWithNames.length/modelsWithIds.length*3), 2),
    text: `${modelsWithNames.length} of ${modelsWithIds.length} models have identifier rules`
  });
}

export const createBucketsSensor = (setInfo: Function, state: any) => {
  const { buckets } = state;
  if (!buckets) return;
  setInfo({ level: !buckets.length ? 0 : 2, text: `${buckets.length} ${plural('bucket',buckets.length)} created` });
};

export const addFilesSensor = (setInfo: Function, state: any) => {
  const { files, bytes } = state;
  console.log({state});

  if (files == null || bytes == null) return;
  setInfo({
    level: files < 10 ? 0 : (bytes < 10000000000 ? 1 : 2),
    text: `${files} ${plural('file',files)}, ${byteFormat(bytes)} stored`
  });
};

export const linkRecordsSensor = (setInfo: Function, state: any) => setInfo(
      { level: 0, text: '200 records created' }
    );

export const createLoadersSensor = (setInfo: Function, state: any) => setInfo(
  { level: 0, text: '1 data loader, last run 2025-02-02' }
);

export const addWorkflowsSensor = (setInfo: Function, state: any) => setInfo(
  { level: 0, text: '1 workflow' }
);

export const runWorkflowsSensor = (setInfo: Function, state: any) => setInfo(
  { level: 0, text: '1 workspace, last run 2025-03-03' }
);
