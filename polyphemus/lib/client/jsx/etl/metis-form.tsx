import React, {
  useState,
  useEffect,
  useCallback,
  useContext,
  useMemo
} from 'react';
import Grid from '@material-ui/core/Grid';
import Tabs from '@material-ui/core/Tabs';
import Tab from '@material-ui/core/Tab';
import Typography from '@material-ui/core/Typography';
import Tooltip from '@material-ui/core/Tooltip';
import Select from '@material-ui/core/Select';
import MenuItem from '@material-ui/core/MenuItem';
import IconButton from '@material-ui/core/IconButton';
import Button from '@material-ui/core/Button';
import AddIcon from '@material-ui/icons/Add';
import DeleteIcon from '@material-ui/icons/Delete';
import CopyIcon from '@material-ui/icons/FileCopy';
import Pagination from '@material-ui/lab/Pagination';
import Switch from '@material-ui/core/Switch';
import FormControlLabel from '@material-ui/core/FormControlLabel';

import {SchemaContext, SchemaProvider} from './schema-context';
import {PickBucket} from 'etna-js/components/metis_exploration';

import {Script, Job} from '../polyphemus';
import AddModel from './add-model';
import { AttributeName, ColumnMap, DefaultItem, FileMatch, FolderPath, TableFormat, useMetisFormStyles, BlankTable } from './metis-form-components';

const SCRIPT_TYPES = [ 'file', 'file_collection', 'data_frame' ];

const SCRIPT_ITEMS = {
  file_match: FileMatch,
  attribute_name: AttributeName,
  folder_path: FolderPath,
  format: TableFormat,
  column_map: ColumnMap,
  blank_table: BlankTable,
  // extracted_columns: ExtractedColumns,
  default: DefaultItem
};

type SchemaType = {
  type?: string;
  enum?: string[]
  const?: string;
}

const defaultFor = (schema_type:SchemaType) => {
  if (schema_type.type == 'string') return '';
  if (schema_type.type == 'array') return [];
  if (schema_type.type == 'object') return {};
  if (schema_type.enum) return schema_type.enum[0];
  if (schema_type.const) return schema_type.const;

  return undefined;
};

const MetisScript = ({
  script,
  num,
  update,
  copy,
  modelName,
  projectName,
  bucketName
}: {
  script: Script;
  num: number;
  update: (newScript: Script | undefined) => void;
  copy: () => void;
  modelName: string;
  projectName: string;
  bucketName: string;
}) => {
  const classes = useMetisFormStyles();
  const {type} = script;

  const handleUpdate = useCallback(
    (newScript: Script | undefined) => {
      update(newScript);
    },
    [update]
  );
  const {schema} = useContext(SchemaContext);


  const defaultScript = useCallback(
    (type: string) => {
      const type_def = schema?.definitions?.[`metis_${type}`];

      if (!type_def) return {};

      const script = type_def.required.reduce(
        (ds:any, req:string) => {
          if (req == 'type') return ds;
          ds[req] = defaultFor(type_def.properties[req]);
          return ds;
        },
        {}
      )

      return script;
    }, [schema]
  );

  const required_fields = schema?.definitions?.[`metis_${type}`]?.required || [];
  const all_fields = Object.keys(schema?.definitions?.[`metis_${type}`]?.properties || {});

  return (
    <Grid className={classes.script} container direction='column'>
      <Typography className={classes.number}>{num}</Typography>
      <Grid
        container
        item
        className={classes.script_header}
        alignItems='center'
      >
        <Grid item className={classes.attribute_name} xs={2}>
          type
        </Grid>
        <Grid item xs={6}>
          <Select displayEmpty value={type || ''} onChange={e => handleUpdate({ type: e.target.value, ...defaultScript(e.target.value as string)})}>
            <MenuItem value='' disabled>
              <em>Select type</em>
            </MenuItem>
            {SCRIPT_TYPES.map((v: string) => (
              <MenuItem key={v} value={v}>
                {v}
              </MenuItem>
            ))}
          </Select>
        </Grid>
        <Grid item xs={4}>
          <Grid container justify='flex-end'>
            <Tooltip title='Copy script'>
              <IconButton onClick={copy}>
                <CopyIcon fontSize='small' />
              </IconButton>
            </Tooltip>
            <Tooltip title='Delete script'>
              <IconButton onClick={() => update(undefined)}>
                <DeleteIcon fontSize='small' />
              </IconButton>
            </Tooltip>
          </Grid>
        </Grid>
      </Grid>
      {
        type && <Grid className={classes.script_items}>
        {
          all_fields.filter( a => a != 'type' ).map(
            (field_name, i) => {
              const value = script[field_name] || '';
              const InputComponent = (SCRIPT_ITEMS[field_name as keyof typeof SCRIPT_ITEMS] || SCRIPT_ITEMS.default) as any;
              return <Grid alignItems='center' key={i} item container>
                <Grid item xs={2}>{field_name}</Grid>
                <Grid item xs={8}>
                  <InputComponent
                    projectName={projectName}
                    bucketName={bucketName}
                    modelName={modelName}
                    type={type}
                    value={value}
                    script={script}
                    update={ (v:any) => handleUpdate({ ...script, [field_name]: v }) }
                    classes={classes}
                  />
                </Grid>
              </Grid>
            }
          )
        }
        </Grid>
      }
    </Grid>
  );
};

const ModelRow = ({
  name,
  children,
  title
}: {
  name: string;
  children: React.ReactNode;
  title?: string;
}) => {
  const classes = useMetisFormStyles();
  return (
    <Grid className={classes.model_row} spacing={1} item container>
      <Tooltip arrow title={title || ''}>
        <Grid
          className={classes.model_row_name}
          item
          container
          justify='flex-end'
          xs={1}
        >
          {name}
        </Grid>
      </Tooltip>
      <Grid item alignItems='center' container xs={11}>
        {children}
      </Grid>
    </Grid>
  );
};

type Model = {
  scripts: Script[];
};

const MetisModel = ({
  config,
  modelName,
  projectName,
  bucketName,
  update
}: {
  config: Model;
  modelName: string;
  projectName: string;
  bucketName: string;
  update: (modelConfig: Model | undefined) => void;
}) => {
  const classes = useMetisFormStyles();
  const {scripts} = config;
  const [pageSize, setPageSize] = useState(5);

  const pages = Math.ceil(scripts.length / pageSize);
  const [page, setPage] = useState(1);

  const page_scripts = scripts.slice((page - 1) * pageSize, page * pageSize);

  const handleUpdateScript = (index: number, newScript: Script | undefined) => {
    const pos = index + (page - 1) * pageSize;
    const newScripts =
      newScript === undefined
        ? scripts.filter((s, j) => j != pos)
        : Object.assign([], scripts, {[pos]: newScript});
    update({...config, scripts: newScripts});
  };

  const handleCopyScript = useCallback(
    (index: number, script: Script) => {
      update({
        ...config,
        scripts: [
          ...scripts.slice(0, index),
          JSON.parse(JSON.stringify(script)),
          ...scripts.slice(index)
        ]
      });
    },
    [config, update, scripts]
  );

  return (
    <Grid className={classes.model} container>
      <ModelRow name='remove'>
        <Button onClick={() => update(undefined)}>Remove model</Button>
      </ModelRow>
      <ModelRow name='scripts'>
        <Button
          onClick={() =>
            update({...config, scripts: [{}, ...scripts]})
          }
        >
          <AddIcon fontSize='small' /> Add Script
        </Button>
        {(pages > 1 || pageSize != 5) && (
          <>
            <Typography className={classes.page_size}>Page size</Typography>
            <Select
              value={pageSize}
              onChange={(e) => setPageSize(e.target.value as number)}
            >
              {[5, 10, 100].map((n) => (
                <MenuItem key={n} value={n}>
                  {n}
                </MenuItem>
              ))}
            </Select>
            <Pagination
              count={pages}
              page={page}
              onChange={(e, v) => setPage(v)}
            />
          </>
        )}
        {page_scripts.map((script: Script | undefined, i: number) => (
          <MetisScript
            key={i}
            script={script}
            num={i + (page - 1) * pageSize + 1}
            modelName={modelName}
            projectName={projectName}
            bucketName={bucketName}
            update={(newScript: Script | undefined) =>
              handleUpdateScript(i, newScript)
            }
            copy={() => handleCopyScript(i + (page - 1) * pageSize, script)}
          />
        ))}
      </ModelRow>
    </Grid>
  );
};

export type Config = {
  bucket_name: string,
  autolink: boolean,
  models: { [modelName: string]: Model };
};

const MetisForm = ({
  config,
  project_name,
  job,
  update
}: {
  project_name: string;
  config: Config;
  job: Job | undefined;
  update: (config: Config) => void;
}) => {
  const classes = useMetisFormStyles();

  const {setSchema} = useContext(SchemaContext);

  useEffect(() => {
    if (job) setSchema(job.schema);
  }, [job]);

  const configModels = config.models || {};

  const modelNames = Object.keys(configModels).sort();

  const bucket_name = config.bucket_name || '';

  const [tab, setTab] = useState(0);

  const modelName = modelNames[tab];
  const [showAddModel, setShowAddModel] = useState(false);

  const handleUpdateModel = useCallback(
    (modelConfig: Model | undefined) => {
      let c = {...config, models: { ...configModels, [modelName]: modelConfig} };
      if (modelConfig === undefined) delete c.models[modelName];

      update(c as Config);
    },
    [config, modelName, update]
  );

  const setBucket = useCallback(
    (bucket_name: string) => {
      let c = {...config, bucket_name };
      update(c as Config);
    },
    [config, update]
  )

  const setAutoLink = useCallback(
    (autolink: boolean) => {
      let c = {...config, autolink };
      update(c as Config);
    },
    [config, update]
  )

  const modelConfig = useMemo(() => {
    return configModels[modelName];
  }, [config, modelName]);

  return (
    <Grid container className={classes.form} direction='column'>
      <Grid item container>
        <Grid item xs={3} style={{paddingLeft:'10px'}}>Autolink Parent Identifiers</Grid>
        <Grid item xs={9}>
          <Switch
            checked={config.autolink}
            onChange={()=>setAutoLink(!config.autolink)}
            color="primary"
          />
        </Grid>
      </Grid>
      <Grid item container>
        <Grid item xs={3} style={{paddingLeft:'10px'}}>Bucket</Grid>
        <Grid item xs={9}>
          <PickBucket setBucket={setBucket} project_name={project_name} bucket={bucket_name} label={null}/>
        </Grid>
      </Grid>
      <Grid item container className={classes.tabs}>
        <Grid item className={classes.tab_scroll}>
          <Tabs
            variant='scrollable'
            indicatorColor='primary'
            scrollButtons='auto'
            value={tab}
            onChange={(e, tab) => setTab(tab)}
          >
            {modelNames.map((modelName) => (
              <Tab label={modelName} key={modelName} />
            ))}
          </Tabs>
        </Grid>
        <Grid item container className={classes.tab_buttons} justify='flex-end'>
          <Tooltip title='Add model'>
            <Button onClick={() => setShowAddModel(true)} startIcon={<AddIcon/>}>
              Add Model
            </Button>
          </Tooltip>
        </Grid>
      </Grid>
      <Grid item className={classes.tab_pane}>
        {modelName != undefined && (
          <MetisModel
            key={modelName}
            modelName={modelName}
            projectName={project_name}
            bucketName={bucket_name}
            config={modelConfig}
            update={handleUpdateModel}
          />
        )}
      </Grid>
      <AddModel
        open={showAddModel}
        close={() => setShowAddModel(false)}
        update={ model_name => update({...config, models: { ...configModels, [model_name]: { scripts: [] } } })}
        excludeModels={Object.keys(config)}
      />
    </Grid>
  );
};

const MetisFormProvider = (props: any) => (
  <SchemaProvider>
    <MetisForm {...props} />
  </SchemaProvider>
);

export default MetisFormProvider;
