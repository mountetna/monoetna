import React, {
  useMemo,
  useCallback,
  useEffect,
  useState,
  useReducer
} from 'react';

import Grid from '@material-ui/core/Grid';
import {makeStyles} from '@material-ui/core/styles';

import AttributeReport from './model_map/attribute_report';
import MapHeading from './model_map/map_heading';
import ModelReport from './model_map/model_report';
import ModelMapGraphic from './model_map/model_map_graphic';

import {selectUser} from 'etna-js/selectors/user-selector';
import {isAdmin} from 'etna-js/utils/janus';
import {useModal} from 'etna-js/components/ModalDialogContainer';
import {requestModels} from 'etna-js/actions/magma_actions';
import {selectTemplate, selectModelNames} from 'etna-js/selectors/magma';
import {fetchProjectsAction} from 'etna-js/actions/janus-actions';
import {selectProjects} from 'etna-js/selectors/janus-selector';
import {projectNameFull} from 'etna-js/utils/janus';
import {useReduxState} from 'etna-js/hooks/useReduxState';
import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import useMagmaActions from './model_map/use_magma_actions';

const mapStyle = makeStyles((theme) => ({
  report: {
    borderLeft: '1px solid #bbb',
    height: 'calc(100vh - 61px)',
    width: 'calc(100vw - 600px)',
    flexWrap: 'nowrap'
  },
  model_map: {
    position: 'relative',
    height: 'calc(100vh - 61px)'
  },
  map: {
    flexBasis: '600px'
  },
  heading: {
    position: 'absolute',
    left: '15px',
    top: '0px',
    height: '48px',
    width: '560px'
  }
}));

const countsReducer = (state, action) => {
  switch (action.type) {
    case 'MODEL_COUNT':
      return {
        ...state,
        [action.model_name]: {
          ...(state[action.model_name] || {}),
          count: action.count
        }
      };
    case 'ATTRIBUTE_COUNT':
      return {
        ...state,
        [action.model_name]: {
          ...(state[action.model_name] || {}),
          attributes: {
            ...state[action.model_name]?.attributes,
            [action.attribute_name]: action.count
          }
        }
      };
    default:
      return state;
  }
};

const ModelMap = ({}) => {
  const [model, setModel] = useState('project');
  const [attribute_name, setAttribute] = useState(null);
  const invoke = useActionInvoker();
  const state = useReduxState();
  const modelNames = useReduxState((state) => selectModelNames(state));
  const template = selectTemplate(state, model);
  const projects = selectProjects(state);
  const user = useReduxState((state) => selectUser(state));

  const attribute =
    template && attribute_name ? template.attributes[attribute_name] : null;
  const updateModel = (model) => {
    setModel(model);
    setAttribute(null);
  };

  const [counts, updateCounts] = useReducer(countsReducer, {});

  useEffect(() => {
    invoke(requestModels());
    invoke(fetchProjectsAction());
  }, []);

  useEffect(() => {
    if (!modelNames.includes(model)) {
      setModel('project');
    }
  }, [modelNames, model]);

  const [width, height] = [600, 600];

  const full_name =
    projectNameFull(projects, CONFIG.project_name) || CONFIG.project_name;

  const classes = mapStyle();

  const isAdminUser = useMemo(() => {
    if (!user || 0 === Object.keys(user).length) return false;

    return isAdmin(user, CONFIG.project_name);
  }, [user, CONFIG.project_name]);

  return (
    <Grid className={classes.model_map} container>
      <Grid item className='map'>
        <MapHeading
          className={classes.heading}
          name='Project'
          title={full_name}
        ></MapHeading>
        <ModelMapGraphic
          width={width}
          height={height}
          selected_models={[model]}
          handler={updateModel}
        />
      </Grid>
      <Grid container direction='column' className={classes.report}>
        <ModelReport
          counts={counts}
          updateCounts={updateCounts}
          key={model}
          model_name={model}
          template={template}
          setAttribute={setAttribute}
          isAdminUser={isAdminUser}
        />
        { attribute && <AttributeReport
            counts={counts[model]}
            attribute={attribute}
            model_name={model}
            isAdminUser={isAdminUser}
            dismiss={() => setAttribute(null)}
          />
        }
      </Grid>
    </Grid>
  );
};

export default ModelMap;
