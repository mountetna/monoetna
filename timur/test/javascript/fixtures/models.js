export const modelsFor = (...model_names) => Object.fromEntries(
  model_names.map(
    model_name => [ model_name, {
      documents: {},
      revisions: {},
      views: {},
      template: require(`./template_${model_name}.json`)
    } ]
  )
);

export const models = modelsFor(
  'prize',
  'monster',
  'labor',
  'project',
  'victim',
  'demographics',
  'wound',
  'aspect',
  'vegetation',
  'habitat'
);
