import React, { useState, useEffect, useCallback, useContext, useMemo } from 'react';
import Grid from '@material-ui/core/Grid';
import Card from '@material-ui/core/Card';
import Tabs from '@material-ui/core/Tabs';
import Tab from '@material-ui/core/Tab';
import Typography from '@material-ui/core/Typography';
import Tooltip from '@material-ui/core/Tooltip';
import TextField from '@material-ui/core/TextField';
import Checkbox from '@material-ui/core/Checkbox';
import CheckBoxOutlineBlankIcon from '@material-ui/icons/CheckBoxOutlineBlank';
import CheckBoxIcon from '@material-ui/icons/CheckBox';
import FormControl from '@material-ui/core/FormControl';
import InputLabel from '@material-ui/core/InputLabel';
import Select from '@material-ui/core/Select';
import MenuItem from '@material-ui/core/MenuItem';
import IconButton from '@material-ui/core/IconButton';
import Button from '@material-ui/core/Button';
import ClearIcon from '@material-ui/icons/Clear';
import AddIcon from '@material-ui/icons/Add';
import DeleteIcon from '@material-ui/icons/Delete';
import CopyIcon from '@material-ui/icons/FileCopy';
import EditIcon from '@material-ui/icons/Edit';
import TuneIcon from '@material-ui/icons/Tune';
import Pagination from '@material-ui/lab/Pagination';

import Dialog from '@material-ui/core/Dialog';
import DialogActions from '@material-ui/core/DialogActions';
import DialogContent from '@material-ui/core/DialogContent';
import DialogContentText from '@material-ui/core/DialogContentText';
import DialogTitle from '@material-ui/core/DialogTitle';

import { ModalProps } from '@material-ui/core/Modal';

import { MagmaContext } from 'etna-js/contexts/magma-context';
import { RedcapContext, RedcapProvider } from './redcap-context';
import {Debouncer} from 'etna-js/utils/debouncer';

import {makeStyles, Theme} from '@material-ui/core/styles';

const useStyles = makeStyles( (theme:Theme) => ({
  form: {
    height: '500px',
    flexWrap: 'nowrap'
  },
  entity: {
    height: '30px',
    width: 'auto',
    background: '#eee',
    padding: '5px',
    borderRadius: '2px 0px 0px 2px'
  },
  filter_entity: {
    '&:not(:last-of-type)': {
      marginRight: '0px',
      borderRadius: '0px'
    },
    height: '30px',
    width: 'auto',
    background: '#fff',
    padding: '5px',
    borderRadius: '0px 2px 2px 0px',
    border: '1px solid #eee',
    marginRight: '5px',
  },
  remove_entity: {
    height: '30px',
    width: 'auto',
    background: '#eee',
    borderRadius: '0px 2px 2px 0px',
    marginRight: '5px'
  },
  tabs: {
    flex: '1 1 auto',
    border: '1px solid #ccc',
    borderBottom: 'none'
  },
  tab_scroll: {
    flexGrow: 0,
    maxWidth: 'calc(100% - 50px)',
    flexBasis: 'calc(100% - 50px)',
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
    "&:not(:last-of-type)": {
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
  attribute_value: {
    width: 'auto',
    flex: '1 1 auto'
  },
  buttons: {
    flex: '0 0 85px'
  },
  field_name: {
    color: theme.palette.secondary.main,
    flex: '1 1 auto',
    marginRight: '5px',
    padding: '0px 10px',
  },
  value: {
    flex: '1 1 auto',
    background: '#fff',
    padding: '0px 5px'
  },
  value_row: {
    borderRadius: '2px',
    border: 'none',
    marginRight: '5px'
  },
  remove_value: {
    width: 'auto',
    flex: '1 1 auto'
  }
}));

const diff = (list1:string[], list2:string[]) => {
  const l2 = new Set(list2);
  return list1.filter( x => !l2.has(x) );
}

const SmallCheckbox = (props:any) => <Checkbox
  icon={<CheckBoxOutlineBlankIcon fontSize="small" />}
  checkedIcon={<CheckBoxIcon fontSize="small" />}
  {...props}/>;

const debounce = (callback:Function, delay:number) => {
  let timer:ReturnType<typeof setTimeout>;
  return (...args:any[]) => {
    clearTimeout(timer);
    timer = setTimeout(() => callback(...args), delay);
  }
}

const ValueRow = ({field_name, value, update, opts}:{
  field_name: string,
  value: any,
  update: (newValue:any) => void,
  opts:any
}) => {
  const classes = useStyles();

  const [localValue, setLocalValue] = useState<typeof value>('');

  useEffect(() => { setLocalValue(value) }, [value]);

  let valueComponent;

  const [debouncer, setDebouncer] = useState(new Debouncer({windowMs: 200}));

  const debouncedUpdate = useCallback(
    (newValue: any) => {
      setLocalValue(newValue);
      debouncer.ready(() => update(newValue));
    },
    [ value, debouncer, update ]
  );

  if (opts === undefined)
    valueComponent = <Typography>{ value }</Typography>;
  else if (opts.type === 'string')
    valueComponent = <TextField fullWidth value={localValue} onChange={(e) => debouncedUpdate(e.target.value)}/>;
  else if (opts.type == 'array')
    valueComponent = <TextField placeholder='Comma-separated list' fullWidth value={localValue.join(', ')} onChange={(e:React.ChangeEvent<HTMLInputElement>) => debouncedUpdate(e.target.value.split(/,\s*/))}/>;
  else if (opts.type == 'boolean')
    valueComponent = <SmallCheckbox checked={value} onChange={(e:React.ChangeEvent<HTMLInputElement>) => update(e.target.checked)}/>;
  else if (opts.enum)
    valueComponent = <Select value={value} onChange={e => update(e.target.value) } >
      {
        opts.enum.map(
           (v:string) => <MenuItem key={v} value={v}>{v}</MenuItem>
        )
      }
    </Select>;

  return <Card className={classes.value_row} variant='outlined' elevation={0}>
    <Grid item container alignItems='center'>
      <Grid item className={classes.field_name}><Typography>{field_name}</Typography></Grid>
      <Grid item className={classes.value}>{ valueComponent }</Grid>
      <Grid item className={classes.remove_value} container justify='flex-end'>
        <Tooltip title='Remove property'>
          <IconButton onClick={ () => update(undefined) }><ClearIcon fontSize='small'/></IconButton>
        </Tooltip>
      </Grid>
    </Grid>
  </Card>
}

const AddProp = ({open,close,update,attribute_value,attribute_props}:{
  open: boolean,
  close: () => void,
  update: (newValue:any) => void,
  attribute_value: { [key: string]: any }|string,
  attribute_props: { [key: string]: any }
}) => {
  const props = diff(
    attribute_props ? Object.keys(attribute_props) : [],
    typeof attribute_value === 'string' ? [ 'redcap_field', 'value' ] : Object.keys(attribute_value)
  );

  const [ newProp, setNewProp ] = useState('');

  return <Dialog open={open} onClose={close}>
    <DialogTitle>Add property</DialogTitle>
    <DialogContent>
      <Select displayEmpty value={newProp} onChange={ (e:React.ChangeEvent<{ value: unknown }>) => setNewProp(e.target.value as string)} >
        <MenuItem value=''><em>None</em></MenuItem>
        {
          props.map(
            prop => <MenuItem key={prop} value={prop}>{prop}</MenuItem>
          )
        }
      </Select>
    </DialogContent>
    <DialogActions>
      <Button disabled={ !newProp } onClick={ () => {
        const opts = attribute_props?.[newProp];
        const newValue = !opts ? undefined :
          opts.type == 'string' ? '' :
          opts.enum ? opts.enum[0] :
          opts.type == 'array' ? [] :
          opts.type == 'boolean' ? true : undefined;

        const oldValue = (typeof attribute_value === 'string') ? { redcap_field: attribute_value, value: 'value' } : attribute_value;

        if (newValue != undefined) update({
          ...oldValue, [newProp]: newValue
        });

        close();
      } } color="secondary">Add</Button>
    </DialogActions>
  </Dialog>
}

const RedcapAttribute = ({att_name, attribute_value, update}:{
  att_name: string,
  attribute_value: any,
  update: (newValue:any) => void
}) => {
  const classes = useStyles();
  const { schema } = useContext(RedcapContext);
  const attribute_props = schema?.definitions?.attribute_value?.properties;

  const [ showAddProp, setShowAddProp ] = useState(false);

  let valueComponent;

  const handleUpdateValue = useCallback((field_name: string, newValue: any) => {
    let v = { ...attribute_value, [field_name]: newValue };
    if (newValue === undefined) delete v[field_name];
    update(v);
  }, [attribute_value, update]);

  if (typeof attribute_value === 'string')
    valueComponent = <TextField placeholder='redcap_field' onChange={ e => update(e.target.value) } value={attribute_value}/>;
  else {
    const sort:{ [key:string]:number|undefined}  = { value: 2, redcap_field: 1 }
    const fieldNames = Object.keys(attribute_value).sort( (a, b) => (sort[b] || 0) - (sort[a] || 0) );

    valueComponent = fieldNames.map(
      field_name => <ValueRow 
        key={field_name}
        field_name={field_name}
        value={attribute_value[field_name]}
        opts={ attribute_props?.[field_name] }
        update={
          (newValue:any) => {
            handleUpdateValue(field_name, newValue);
          }
        }
      />
    )
  }

  return <Grid key={att_name} item container alignItems='center'>
    <Tooltip placement='left' title={att_name}>
      <Grid className={ classes.attribute_name} item xs={2}>
        {att_name}
      </Grid>
    </Tooltip>
    <Grid item container className={classes.attribute_value}>
      {valueComponent}
    </Grid>
    <Grid item>
      <Tooltip title='Add property'><IconButton onClick={() => setShowAddProp(true)}><EditIcon fontSize='small'/></IconButton></Tooltip>
      <Tooltip title='Remove attribute'><IconButton onClick={() => update(undefined)}><ClearIcon fontSize='small'/></IconButton></Tooltip>
    </Grid>
    <AddProp
      attribute_value={attribute_value}
      open={showAddProp}
      close={() => setShowAddProp(false)}
      update={update}
      attribute_props={attribute_props}
    />
  </Grid>
}

const AddAttribute = ({open,close,update,script,modelName}:{
  open:boolean,
  close:()=>void,
  update:Function,
  script:Script,
  modelName:string
}) => {
  const [ newAttribute, setNewAttribute ] = useState('');

  const { models } = useContext(MagmaContext);

  const attribute_names = Object.keys(models[modelName]?.template?.attributes || {});

  return <Dialog open={open} onClose={close}>
    <DialogTitle>Add Attribute</DialogTitle>
    <DialogContent>
      <Select displayEmpty value={newAttribute} onChange={ e => setNewAttribute(e.target.value as string)} >
        <MenuItem value=''><em>None</em></MenuItem>
        {
          attribute_names.map( att_name => <MenuItem key={att_name} value={att_name}>{att_name}</MenuItem>)
        }
      </Select>
    </DialogContent>
    <DialogActions>
      <Button disabled={ !newAttribute } onClick={ () => {
        update({ ...script, attributes: { ...script.attributes, [newAttribute]: '' }});

        close();
      } } color="secondary">Add</Button>
    </DialogActions>
  </Dialog>
}

type Entity = string|{[key:string]:string};

const isFilteredEntity = (e:Entity) => typeof e !== 'string';

const entityName = (e:Entity):string => typeof e === 'string' ? e : Object.keys(e)[0];
const entityValue = (e:Entity):string => Object.values(e)[0];

const AddEntity = ({open,close,update,each}:{
  open:boolean,
  close:()=>void,
  update: (newEach:Entity[]|undefined) => void,
  each:Entity[]|undefined
}) => {
  const [ newEntity, setNewEntity ] = useState('');

  each = each || [];

  const { schema } = useContext(RedcapContext);
  const entity_names = diff(
    Object.keys(schema?.definitions?.each_entity?.properties || {}),
    each.map(entityName)
  );

  return <Dialog open={open} onClose={close}>
    <DialogTitle>Add Entity</DialogTitle>
    <DialogContent>
      <Select displayEmpty value={newEntity} onChange={ e => setNewEntity(e.target.value as string)} >
        <MenuItem value=''><em>None</em></MenuItem>
        {
          entity_names.map( entity_name => <MenuItem key={entity_name} value={entity_name}>{entity_name}</MenuItem>)
        }
      </Select>
    </DialogContent>
    <DialogActions>
      <Button disabled={ !newEntity } onClick={ () => {
        let newEach = each !== undefined ? [ ...each, newEntity ] : [ newEntity ];

        update(newEach);

        close();
      } } color="secondary">Add</Button>
    </DialogActions>
  </Dialog>
}
const RedcapEntity = ({each, allowNull=false, update}:{
  each: Entity[]|undefined,
  allowNull?: boolean,
  update: (newEach:Entity[]|undefined) => void
}) => {
  const classes = useStyles();

  const updateItem = useCallback(
    (i, entity) => update(Object.assign([], each, {[i]: entity})), [update, each]
  );
  
  const [ showAddEntity, setShowAddEntity ] = useState(false);

  return <Grid item container alignItems='center'>
    {
      !each ? <em>inherited</em> :
      each.map( (entity,i) => <Grid item key={i}>
        <Grid container alignContent='center'>
          <Grid item container alignContent='center' className={classes.entity}>{ entityName(entity) }</Grid>
          <Grid item container alignContent='center' className={classes.filter_entity}>
            <Tooltip title='Filter entity'><IconButton
              color={ isFilteredEntity(entity) ? 'secondary' : 'default' }
              onClick={ () => isFilteredEntity(entity) ? updateItem(i, entityName(entity)) : updateItem(i, { [entity as string]: '' }) }
              size='small'><TuneIcon fontSize='small'/></IconButton></Tooltip>
            {
              isFilteredEntity(entity) ? 
                <TextField
                  onChange={ e => updateItem(i, { [ entityName(entity) ]: e.target.value }) }
                  placeholder='Filter'
                  value={entityValue(entity)}/> : null
            }
          </Grid>
          {
            (entityName(entity) != 'record' || each.length == 1 && allowNull) &&
            <Grid item container alignContent='center' className={classes.remove_entity}>
              <Tooltip title='Remove entity'><IconButton onClick={
                () => each.length == 1 ? update(undefined) : update(each.filter( (e,j) => j != i ))
                } size='small'><ClearIcon fontSize='small'/></IconButton></Tooltip>
            </Grid>
          }
        </Grid>
      </Grid>)
    }
    {
      (!each || each.length < 4) &&
      <Tooltip title='Add entity'><IconButton onClick={
        ()=> {
          if (each === undefined || each.length == 0) update(['record']);
          else setShowAddEntity(true);
        }
      } size='small'><AddIcon fontSize='small'/></IconButton></Tooltip>
    }
    <AddEntity
      open={showAddEntity}
      close={() => setShowAddEntity(false)}
      update={update}
      each={each}
    />
  </Grid>
}

type Script = any;

const RedcapScript = ({script, num, update, copy, modelName}:{
  script: Script,
  num: number,
  update: (newScript:Script|undefined) => void,
  copy: () => void,
  modelName: string
}) => {
  const classes = useStyles();
  const { attributes, each } = script;

  const [ showAddAttribute, setShowAddAttribute ] = useState(false);

  const handleUpdate = useCallback((newScript: Script | undefined) => {
    update(newScript);
  }, [update]);

  const updateAttributes = useCallback((att_name, newValue) => {
    let s = { ...script, attributes: { ...attributes, [att_name]: newValue } }
    if (newValue === undefined) delete s.attributes[att_name];
    handleUpdate(s);
  }, [ script, attributes, handleUpdate ]);

  const handleUpdateRedcapEntity = useCallback((newEach: Entity[]|undefined) => {
    const s = { ...script, each: newEach };
    if (newEach === undefined) delete s.each;
    handleUpdate(s);
  }, [handleUpdate, each, script]);

  return <Grid className={ classes.script } container direction='column'>
    <Typography className={classes.number}>{num}</Typography>
    <Grid container item className={ classes.script_header} alignItems='center' >
      <Grid item className={ classes.attribute_name } xs={2}>each</Grid>
      <Grid item xs={6}>
        <RedcapEntity each={each} allowNull={true} update={handleUpdateRedcapEntity} />
      </Grid>
      <Grid item xs={4}>
        <Grid container justify='flex-end'>
          <Tooltip title='Add attribute'><IconButton onClick={ () => setShowAddAttribute(true) }><AddIcon fontSize='small'/></IconButton></Tooltip>
          <Tooltip title='Copy script'><IconButton onClick={ copy }><CopyIcon fontSize='small'/></IconButton></Tooltip>
          <Tooltip title='Delete script'><IconButton onClick={ () => update(undefined) }><DeleteIcon fontSize='small'/></IconButton></Tooltip>
        </Grid>
      </Grid>
    </Grid>
    {
      Object.keys(attributes).map(
        att_name => <RedcapAttribute
          key={att_name}
          att_name={att_name}
          attribute_value={attributes[att_name]}
          update={ (newValue:any) => updateAttributes(att_name, newValue) } />
      )
    }
    <AddAttribute
      modelName={modelName}
      script={script}
      open={showAddAttribute}
      close={() => setShowAddAttribute(false)}
      update={update}
    />
  </Grid>
}

const ModelRow = ({name,children}:{
  name: string,
  children: React.ReactNode
}) => {
  const classes = useStyles();
  return <Grid className={classes.model_row} spacing={1} item container>
    <Grid className={classes.model_row_name} item container justify='flex-end' xs={1} >{name}</Grid>
    <Grid item alignItems='center' container xs={11}>{ children }</Grid>
  </Grid>
}

type Model = {
  scripts: Script[],
  each: Entity[],
  invert?: boolean,
  identifier_fields?: string[]
};

const RedcapModel = ({config,modelName, update}:{
  config: Model,
  modelName: string,
  update: (modelConfig:Model|undefined) => void
}) => {
  const classes = useStyles();
  const { each, invert, scripts, identifier_fields } = config;
  const [ pageSize, setPageSize ] = useState(5);

  const pages = Math.ceil(scripts.length / pageSize);
  const [ page, setPage ] = useState(1);

  const page_scripts = scripts.slice((page-1)*pageSize, page*pageSize);

  const handleUpdateScript = (index: number, newScript: Script|undefined) => {
    const pos = index + (page-1)*pageSize;
    const newScripts = (newScript === undefined) ?
      scripts.filter((s,j) => j != pos) :
      Object.assign([], scripts, {[pos]: newScript});
    update({ ...config, scripts: newScripts });
  }

  const handleCopyScript = useCallback((index: number, script: Script) => {
    update({
      ...config,
      scripts: [ ...scripts.slice(0,index),
        JSON.parse(JSON.stringify(script)),
        ...scripts.slice(index)
      ]
    })
  }, [config, update, scripts])

  return <Grid className={classes.model} container>
    <ModelRow name='remove'><Button onClick={() => update(undefined) }>Remove model</Button></ModelRow>
    <ModelRow name='each'><RedcapEntity each={each} update={ (newEach:Entity[]|undefined) => {
      const c = {...config, each: newEach };
      if (newEach === undefined) delete c.each;
      update(c as Model);
    } }/></ModelRow>
    <ModelRow name='invert'><SmallCheckbox size='small' checked={ !!invert } onChange={ (e:React.ChangeEvent<HTMLInputElement>) => update({ ...config, invert: e.target.checked }) }/></ModelRow>
    <ModelRow name='identifier_fields'><TextField placeholder='Comma-separated list' value={ (identifier_fields||[]).join(', ') } onChange={ e => update({ ...config, identifier_fields: e.target.value.split(/,\s*/) }) } /></ModelRow>
    <ModelRow name='scripts'>
      <Button onClick={ () => update({ ...config, scripts: [ { attributes: {} }, ...scripts ]}) }><AddIcon fontSize='small'/> Add Script</Button>
      { (pages > 1 || pageSize != 5) && <>
          <Typography className={classes.page_size}>Page size</Typography>
          <Select value={pageSize} onChange={e => setPageSize(e.target.value as number)}>
            {
              [ 5, 10, 100 ].map( n => <MenuItem key={n} value={n}>{n}</MenuItem>)
            }
          </Select>
          <Pagination count={ pages } page={page} onChange={ (e,v) => setPage(v) }/>
        </>
      }
        { page_scripts.map(
          (script: Script | undefined, i: number) => 
          <RedcapScript
            key={i}
            script={script}
            num={i+(page-1)*pageSize+1}
            modelName={modelName}
            update={
              (newScript:Script|undefined) => 
                handleUpdateScript(i, newScript)
            }
            copy={() => handleCopyScript(i, script)}
          />
      )}
    </ModelRow>
  </Grid>
};

const AddModel = ({open,close,update,config}:{
  open: boolean,
  close: ()=>void,
  update: (config:Config) => void,
  config: Config
}) => {
  const [ newModel, setNewModel ] = useState('');

  const { models } = useContext(MagmaContext);

  const model_names = diff( Object.keys(models), Object.keys(config) )

  return <Dialog open={open} onClose={close}>
    <DialogTitle>Add Model</DialogTitle>
    <DialogContent>
      <Select displayEmpty value={newModel} onChange={ e => setNewModel(e.target.value as string)} >
        <MenuItem value=''><em>None</em></MenuItem>
        {
          model_names.map( att_name => <MenuItem key={att_name} value={att_name}>{att_name}</MenuItem>)
        }
      </Select>
    </DialogContent>
    <DialogActions>
      <Button disabled={ !newModel } onClick={ () => {
        update({ ...config, [newModel]: { each: ['record'], scripts: [] } })

        close();
      } } color="secondary">Add</Button>
    </DialogActions>
  </Dialog>
}

export type Config = {
  [modelName: string]: Model
}

const RedcapForm = ({config, project_name, job, update}:{
  project_name: string,
  config: Config,
  job: Job|undefined,
  update: (config:Config) => void
}) => {
  const classes = useStyles();

  const { setSchema } = useContext(RedcapContext);

  useEffect( () => {
    if (job) setSchema(job.schema);
  }, [ job ])

  const modelNames = Object.keys(config);

  const [ tab, setTab ] = useState(0);

  const modelName = modelNames[tab];
  const [ showAddModel, setShowAddModel ] = useState(false);

  const handleUpdateModel = useCallback((modelConfig: Model|undefined) => {
    let c = { ...config, [modelName]: modelConfig};
    if (modelConfig === undefined) delete c[modelName];
    update(c as Config);
  }, [config, modelName, update]);

  const modelConfig = useMemo(() => {
    return config[modelName];
  }, [config, modelName]);

  return <Grid container className={classes.form} direction='column'>
    <Grid item container className={classes.tabs}>
      <Grid item className={classes.tab_scroll}>
        <Tabs
          variant="scrollable"
          indicatorColor="primary"
          scrollButtons="auto"
          value={tab}
          onChange={(e,tab) => setTab(tab)}
        >
        {
          modelNames.map( modelName => <Tab label={modelName} key={modelName}/>)
        }
        </Tabs>
      </Grid>
      <Grid item container className={classes.tab_buttons} justify='flex-end'>
        <Tooltip title='Add model'><IconButton onClick={ () => setShowAddModel(true) }><AddIcon/></IconButton></Tooltip>
      </Grid>
    </Grid>
    <Grid item className={classes.tab_pane}>
      {
        modelName != undefined && <RedcapModel
          key={modelName}
          modelName={modelName}
          config={modelConfig}
          update={handleUpdateModel}/>
      }
    </Grid>
    <AddModel
      open={showAddModel}
      close={() => setShowAddModel(false)}
      update={update}
      config={config}
    />
  </Grid>
}

const RedcapFormProvider = (props:any) => <RedcapProvider>
  <RedcapForm {...props}/>
</RedcapProvider>

export default RedcapFormProvider;
