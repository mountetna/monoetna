import React, {useState, useContext, useCallback, useEffect, useMemo} from 'react';

import ModelMapGraphic from '../model_map/model_map_graphic';
import { Model, Models } from 'etna-js/models/magma-model';
import { MagmaContext } from 'etna-js/contexts/magma-context';
import {plural} from 'etna-js/utils/format';
import {requestModels} from 'etna-js/actions/magma_actions';
import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import MagmaLink from '../magma_link';
import {packParams, unpackParams} from '../../utils/query_uri_params';
import {columnsForModel} from '../../contexts/query/query_column_context';

import { useTheme, makeStyles } from '@material-ui/core/styles';
import Grid from '@material-ui/core/Grid';
import Popper from '@material-ui/core/Popper';
import Card from '@material-ui/core/Card';
import CardContent from '@material-ui/core/CardContent';
import Tooltip from '@material-ui/core/Tooltip';
import IconButton from '@material-ui/core/IconButton';
import ChevronRightIcon from '@material-ui/icons/ChevronRight';
import ExpandMoreIcon from '@material-ui/icons/ExpandMore';
import YoutubeSearchedForIcon from '@material-ui/icons/YoutubeSearchedFor';

const useStyles = makeStyles((theme) => ({
  model_view: {
    width: '200px',
    height: '26px',
    fontFamily: 'sans-serif',
    marginTop: '5px',
    background: 'white',
    fontSize: '0.75rem',
    border: '1px solid #0c0',
    borderRadius: '1px',
    boxShadow: '0 0 4px 0 #bbb',
    display: 'flex',
    alignItems: 'center'
  },
  id_list: {
    width: '200px',
    height: '104px',
    overflowY: 'auto',
    display: 'flex',
    flexDirection: 'column',
    border: '1px solid #ccc'
  },
  id: {
    display: 'flex',
    alignItems: 'center',
    width: 'calc(100% - 4px)',
    padding: '2px',
    fontFamily: 'monospace',
    borderBottom: '1px solid #eee'
  },
  model_name: {
    flex: '1 1 auto',
    textAlign: 'right',
    paddingRight: '5px'
  }
}));

const ModelMiniReport = ({modelEl, modelName, open}:{
  modelName: string,
  modelEl: HTMLElement | null,
  open: boolean
}) => {
  const [fold, setFold] = useState(true);
  const classes = useStyles();

  useEffect( () => {
    setFold(true);
  }, [modelName] );

  const { counts, models: modelsObj } = useContext(MagmaContext);

  const models = useMemo( () => new Models(modelsObj), [ modelsObj ] );

  const model = models.model(modelName);

  const setQueryLocation = useCallback(
    async () => {
      const columns = columnsForModel(model as Model, models);
      const paramString = await packParams({
        rootModel: modelName,
        recordFilters: [],
        orRecordFilterIndices: [],
        columns
      });
      window.open(`/${CONFIG.project_name}/query#${paramString}`, '_blank');
    }, [models, model, modelName] );

  return <Popper anchorEl={modelEl} open={open}>
    <Card elevation={0}>
      <Grid className={ classes.model_view }>
        <Grid className={classes.model_name}>
          { modelName ? `${counts[modelName] == undefined ? '??' : counts[modelName]} ${plural('record',counts[modelName])}` : null }
        </Grid>
        <Tooltip title={ model?.isTable ? `Browse via ${model.parent} model` : counts[modelName] == 0 ? 'No records' : 'Show ids' }>
          <span>
          <IconButton size='small' disabled={ model?.isTable || counts[modelName] == 0 } onClick={ () => setFold(!fold) }>
            { fold ? <ChevronRightIcon fontSize='small'/> : <ExpandMoreIcon fontSize='small'/> }
          </IconButton>
          </span>
        </Tooltip>
        <Tooltip title={ model && counts[modelName] == 0 ? 'No records' : 'Query and Download' }>
          <span>
            <IconButton size='small' disabled={ counts[modelName] == 0 } onClick={
              setQueryLocation
            }>
              <YoutubeSearchedForIcon/>
            </IconButton>
          </span>
        </Tooltip>
      </Grid>
      {
        !fold && open && modelName && <Grid className={ classes.id_list }>
        {
          Object.keys(modelsObj[modelName].documents).map(
            id => <Grid key={id} className={classes.id}>
              <MagmaLink link={id} model={modelName} />
            </Grid>
          )
        }
        </Grid>
      }
    </Card>
  </Popper>
}

type ModelState = [ modelName: string|null, modelEl: HTMLElement|null ];

const DashboardMap = ({modelState, setModelState}:{
  modelState: ModelState,
  setModelState: (s: ModelState) => void
}) => {
  const [modelName, modelEl] = modelState;
  const invoke = useActionInvoker();

  useEffect(() => {
    invoke(requestModels());
  }, []);

  const updateModel = useCallback( (newModelName, newModelEl) => {
    if (modelName && modelName == newModelName) {
      setModelState([null, null]);
    } else {
      setModelState([newModelName, newModelEl]);
    }
  }, [modelName]);

  const showPopper: boolean = Boolean(modelEl);

  return <Grid>
    <ModelMapGraphic
      width={600}
      height={600}
      handler={updateModel}
      selected_models={modelName ? [modelName] : []}
      disabled_models={[]}
    />
    <ModelMiniReport
      open={showPopper}
      modelEl={modelEl}
      modelName={modelName as string}
      />
  </Grid>
}

export default DashboardMap;
