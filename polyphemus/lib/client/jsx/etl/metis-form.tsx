import React, {
  useState,
  useEffect,
  useCallback,
  useContext,
  useMemo
} from 'react';
import Grid from '@material-ui/core/Grid';
import Card from '@material-ui/core/Card';
import Tabs from '@material-ui/core/Tabs';
import Tab from '@material-ui/core/Tab';
import Typography from '@material-ui/core/Typography';
import Tooltip from '@material-ui/core/Tooltip';
import TextField from '@material-ui/core/TextField';
import Select from '@material-ui/core/Select';
import MenuItem from '@material-ui/core/MenuItem';
import IconButton from '@material-ui/core/IconButton';
import Button from '@material-ui/core/Button';
import AddIcon from '@material-ui/icons/Add';
import DeleteIcon from '@material-ui/icons/Delete';
import CopyIcon from '@material-ui/icons/FileCopy';
import Pagination from '@material-ui/lab/Pagination';

import {SchemaContext, SchemaProvider} from './schema-context';
import AddDialog from 'etna-js/components/add-dialog';

import {makeStyles, Theme} from '@material-ui/core/styles';

import {Job} from '../polyphemus';
import AddModel from './add-model';
import {diff} from '../utils/list';
import SelectAttribute from '../select-attribute';

const useStyles = makeStyles((theme: Theme) => ({
  form: {
    height: '500px',
    flexWrap: 'nowrap'
  },
  tabs: {
    flex: '1 1 auto',
    border: '1px solid #ccc',
    borderBottom: 'none'
  },
  tab_scroll: {
    flexGrow: 0,
    maxWidth: 'calc(100% - 50px)',
    flexBasis: 'calc(100% - 50px)'
  },
  tab_buttons: {
    flex: '1 1 auto',
    width: 'auto'
  },
  tab_pane: {
    flex: '1 1 auto',
    border: '1px solid #ccc',
    minHeight: '500px',
    overflow: 'hidden',
    overflowY: 'auto'
  },
  model: {
    padding: '5px'
  },
  model_row: {
    minHeight: '47px',
    '&:not(:last-of-type)': {
      borderBottom: '1px solid #eee',
      marginBottom: '5px'
    }
  },
  model_row_name: {
    marginTop: '7px'
  },
  page_size: {
    color: '#666',
    padding: '0px 15px'
  },
  script: {
    background: '#eee',
    marginTop: '15px',
    border: '1px solid #ccc',
    borderRadius: '2px',
    position: 'relative'
  },
  script_header: {
    borderBottom: '1px solid #ccc',
    background: '#ccc'
  },
  number: {
    color: '#888',
    position: 'absolute',
    transform: 'translate(-100%, 0)',
    left: '-10px',
    top: '10px'
  },
  attribute_name: {
    overflow: 'hidden',
    textOverflow: 'ellipsis',
    padding: '0px 5px'
  },
  script_items: {
    padding: '5px'
  },
  strikeout: {
    '&:hover': {
      textDecoration: 'line-through',
      cursor: 'pointer',
      background: '#eee'
    }
  }
}));

const SCRIPT_TYPES = [ 'file', 'file_collection', 'table' ];

type ScriptItem = {
  type: string;
  projectName: string;
  modelName: string;
  value: string;
  update: Function;
};

const SCRIPT_ITEMS = {
  bucket_name: ({value,update,projectName}) => <TextField
    placeholder={`Metis bucket_name for ${projectName}`}
    fullWidth
    value={value}
    onChange={
      (e: React.ChangeEvent<HTMLInputElement>) => update(e.target.value.split(/,\s*/))
    }
  />,
  folder_path: ({value,update}) => <TextField
    placeholder='Metis folder path'
    fullWidth
    value={value}
    onChange={
      (e: React.ChangeEvent<HTMLInputElement>) => update(e.target.value.split(/,\s*/))
    }
  />,
  file_match: ({value,update}) => <TextField
    placeholder='Regular expression matching file'
    fullWidth
    value={value}
    onChange={
      (e: React.ChangeEvent<HTMLInputElement>) => update(e.target.value.split(/,\s*/))
    }
  />,
  attribute_name: ({type,value,update,modelName}:ScriptItem) => <SelectAttribute
    value={value}
    update={update}
    modelName={modelName}
    filter={ a => a.attribute_type == type }
  />,
  folder_path: ({value,update}) => <TextField
    placeholder='Metis folder path'
    fullWidth
    value={value}
    onChange={
      (e: React.ChangeEvent<HTMLInputElement>) => update(e.target.value.split(/,\s*/))
    }
  />,
  format: ({value,update}) => <Select displayEmpty value={value} onChange={e => update(e.target.value)}>
      <MenuItem value=''><em>Table format</em></MenuItem>
      <MenuItem value='tsv'>tsv</MenuItem>
      <MenuItem value='csv'>csv</MenuItem>
  </Select>,
  column_map: ({value,update}) => {
    const classes = useStyles();
    return <Grid>
      <AddDialog
        title='Map a column'
        content={
          <>Map a column name in the table to another name (e.g., an attribute_name)</>
        }
        buttonText='ADD MAPPING'
        update={ (column_name, mapped_column_name) => update(
          { ...value, [column_name]: mapped_column_name }
        )}
        mask={ v => v }
        mask2={ v => v }
        placeholders={[ 'Column Name', 'Mapped Column Name' ]}
      />
      <Grid item container>
        {
          Object.keys(value).map(
            column_name => <Grid
              key={column_name} item container
              className={classes.strikeout}
              onClick={
                () => {
                  let { [column_name]:_, ...newValue } = value;
                  update(newValue);
                }
              }>
              <Grid item xs={2}>{column_name}</Grid>
              <Grid item xs={8}>{value[column_name]}</Grid>
            </Grid>
          )
        }
      </Grid>
    </Grid>
  },
  extracted_columns: ({value,update}) => <TextField
    placeholder='Comma-separated list'
    fullWidth
    value={Array.isArray(value) ? value.join(', ') : value}
    onChange={
      (e: React.ChangeEvent<HTMLInputElement>) => update(e.target.value.split(/,\s*/))
    }
  />,
  default: ({value, update}:ScriptItem) => <TextField
    fullWidth
    value={value}
    onChange={ e => update(e.target.value as string) }
  />
};

const defaultFor = schema_type => {
  if (schema_type.type == "string") return '';
  if (schema_type.type == "array") return [];
  if (schema_type.type == "object") return {};
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
  projectName
}: {
  script: Script;
  num: number;
  update: (newScript: Script | undefined) => void;
  copy: () => void;
  modelName: string;
  projectName: string;
}) => {
  const classes = useStyles();
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

      const d = type_def.required.reduce(
        (ds, req) => {
          if (req == 'type') return ds;
          ds[req] = defaultFor(type_def.properties[req]);
          return ds;
        },
        {}
      )

      console.log({d, type_def});

      return d;
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
          <Select displayEmpty value={type || ''} onChange={e => handleUpdate({ type: e.target.value, ...defaultScript(e.target.value)})}>
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
              const InputComponent = SCRIPT_ITEMS[field_name] || SCRIPT_ITEMS.default;
              return <Grid alignItems='center' key={i} item container>
                <Grid item xs={2}>{field_name}</Grid>
                <Grid item xs={8}>
                  <InputComponent
                    projectName={projectName}
                    modelName={modelName}
                    type={type}
                    value={value}
                    update={ v => handleUpdate({ ...script, [field_name]: v }) }
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
  const classes = useStyles();
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
  update
}: {
  config: Model;
  modelName: string;
  projectName: string;
  update: (modelConfig: Model | undefined) => void;
}) => {
  const classes = useStyles();
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
  [modelName: string]: Model;
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
  const classes = useStyles();

  const {setSchema} = useContext(SchemaContext);

  useEffect(() => {
    if (job) setSchema(job.schema);
  }, [job]);

  const modelNames = Object.keys(config).sort();

  const [tab, setTab] = useState(0);

  const modelName = modelNames[tab];
  const [showAddModel, setShowAddModel] = useState(false);

  const handleUpdateModel = useCallback(
    (modelConfig: Model | undefined) => {
      let c = {...config, [modelName]: modelConfig};
      if (modelConfig === undefined) delete c[modelName];
      update(c as Config);
    },
    [config, modelName, update]
  );

  const modelConfig = useMemo(() => {
    return config[modelName];
  }, [config, modelName]);

  return (
    <Grid container className={classes.form} direction='column'>
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
            <IconButton onClick={() => setShowAddModel(true)}>
              <AddIcon />
            </IconButton>
          </Tooltip>
        </Grid>
      </Grid>
      <Grid item className={classes.tab_pane}>
        {modelName != undefined && (
          <MetisModel
            key={modelName}
            modelName={modelName}
            projectName={project_name}
            config={modelConfig}
            update={handleUpdateModel}
          />
        )}
      </Grid>
      <AddModel
        open={showAddModel}
        close={() => setShowAddModel(false)}
        update={ model_name => update({...config, [model_name]: { scripts: [] } })}
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
