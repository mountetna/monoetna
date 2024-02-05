import {useDispatch} from 'react-redux';
import React, {useState, useCallback, useMemo, useEffect} from 'react';
import SelectProjectModelDialog from '../select_project_model';
import {requestAnswer} from 'etna-js/actions/magma_actions';
import {getDocuments} from 'etna-js/api/magma_api';
import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';

import Grid from '@material-ui/core/Grid';
import Typography from '@material-ui/core/Typography';
import {makeStyles} from '@material-ui/core/styles';
import Menu from '@material-ui/core/Menu';
import MenuItem from '@material-ui/core/MenuItem';
import MapHeading from './map_heading';
import Button from '@material-ui/core/Button';
import IconButton from '@material-ui/core/IconButton';
import MoreVertIcon from '@material-ui/icons/MoreVert';
import Tooltip from '@material-ui/core/Tooltip';
import LinearProgress from '@material-ui/core/LinearProgress';
import Chip from '@material-ui/core/Chip';
import AddIcon from '@material-ui/icons/Add';
import LinkIcon from '@material-ui/icons/Link';
import DeleteIcon from '@material-ui/icons/Delete';
import SwapHorizIcon from '@material-ui/icons/SwapHoriz';
import LibraryAddIcon from '@material-ui/icons/LibraryAdd';
import FileCopyIcon from '@material-ui/icons/FileCopy';

import AddAttributeModal from './add_attribute_modal';
import AddLinkModal from './add_link_modal';
import ReparentModelModal from './reparent_model_modal';
import AddModelModal from './add_model_modal';
import CopyModelModal from './copy_model_modal';
import RemoveModelModal from './remove_model_modal';
import {useModal} from 'etna-js/components/ModalDialogContainer';
import {selectModels} from 'etna-js/selectors/magma';
import {useReduxState} from 'etna-js/hooks/useReduxState';
import {
  addAttributes,
  addLink,
  removeModel,
  reparentModel,
  addModel
} from '../../api/magma_api';
import {
  addTemplatesAndDocuments,
  removeModelAction
} from 'etna-js/actions/magma_actions';
import {isLeafModel} from '../../utils/edit_map';
import useMagmaActions from './use_magma_actions';
import ModelAttributesTable from './model_attributes_table';

const actionStyles = makeStyles((theme) => ({
  addBtn: {
    margin: '0.5rem 1.0rem 0.5rem 0rem'
  }
}));

const ManageModelActions = ({
  handleAddAttribute,
  handleAddLink,
  handleAddModel,
  handleCopyModel,
  handleRemoveModel,
  handleReparentModel,
  isLeaf,
  canReparent,
  determiningCanReparent,
  modelName
}) => {
  const classes = actionStyles();
  const [ showCopyModelModal, setShowCopyModelModal ] = useState(false);
  const [ showAddModelModal, setShowAddModelModal ] = useState(false);
  const [ showAddAttributeModal, setShowAddAttributeModal ] = useState(false);
  const [ showAddLinkModal, setShowAddLinkModal ] = useState(false);
  const [ showReparentModelModal, setShowReparentModelModal ] = useState(false);
  const [ showRemoveModelModal, setShowRemoveModelModal ] = useState(false);

  const reparentTooltip = determiningCanReparent ? 'Determining if reparenting is possible' :
    canReparent ? 'Reparent Model' : 'Cannot reparent a model containing records'

  return (
    <Grid>
      <Tooltip title='Add Attribute' aria-label='Add Attribute'>
        <Button
          className={classes.addBtn}
          startIcon={<AddIcon />}
          onClick={() => setShowAddAttributeModal(true) }
        >
          Attribute
        </Button>
      </Tooltip>
      <AddAttributeModal
        onClose={ () => setShowAddAttributeModal(false) }
        open={showAddAttributeModal}
        onSave={handleAddAttribute} />
      <Tooltip title='Add Link' aria-label='Add Link'>
        <Button
          className={classes.addBtn}
          startIcon={<LinkIcon />}
          onClick={() => setShowAddLinkModal(true)}
        >
          Link
        </Button>
      </Tooltip>
      <AddLinkModal
        onClose={ () => setShowAddLinkModal(false) }
        open={showAddLinkModal}
        onSave={handleAddLink} />
      <Tooltip title='Add Model' aria-label='Add Model'>
        <Button
          className={classes.addBtn}
          startIcon={<LibraryAddIcon />}
          onClick={() => setShowAddModelModal(true)}
        >
          Model
        </Button>
      </Tooltip>
      <AddModelModal
        modelName={modelName}
        onClose={ () => setShowAddModelModal(false) }
        open={showAddModelModal}
        onSave={handleAddModel} />
      <Tooltip title='Copy Model Attributes' aria-label='Copy Model Attributes'>
        <Button
          className={classes.addBtn}
          startIcon={<FileCopyIcon />}
          onClick={() => setShowCopyModelModal(true)}
        >
          Copy Attributes
        </Button>
      </Tooltip>
      <CopyModelModal
        modelName={modelName}
        onClose={ () => setShowCopyModelModal(false) }
        open={showCopyModelModal}
        onSave={handleCopyModel} />
      <Tooltip title={reparentTooltip} aria-label={reparentTooltip}>
        <span>
          <Button
            className={classes.addBtn}
            startIcon={<SwapHorizIcon />}
            disabled={!canReparent}
            onClick={() => setShowReparentModelModal(true)}
          >
            Reparent Model
          </Button>
        </span>
      </Tooltip>
      <ReparentModelModal
        modelName={modelName}
        onClose={ () => setShowReparentModelModal(false) }
        open={showReparentModelModal}
        onSave={handleReparentModel}
      />
      {isLeaf && <>
        <Tooltip title='Remove Model' aria-label='Remove Model'>
          <Button
            className={classes.addBtn}
            startIcon={<DeleteIcon />}
            onClick={() => setShowRemoveModelModal(true)}
          >
            Remove Model
          </Button>
        </Tooltip>
        <RemoveModelModal
          modelName={modelName}
          onClose={ () => setShowRemoveModelModal(false) }
          open={showRemoveModelModal}
          onSave={handleRemoveModel}
        />
      </>}
    </Grid>
  );
};


const reportStyles = makeStyles((theme) => ({
  model_report: {
    flex: '1 1 auto',
    display: 'flex',
    flexDirection: 'column',
    padding: '0px 10px',
    overflow: 'scroll'
  },
  attributes: {
    flex: '1 1 auto'
  },
  heading: {
    marginBottom: '10px'
  },
  report_row: {
    borderBottom: '1px solid #eee'
  },
  name: {
    color: 'gray'
  },
  title: {
    color: 'forestgreen'
  },
  type: {
    width: '25%',
    paddingRight: '10px'
  },
  count: {
    width: '20%'
  }
}));

const ModelReport = ({
  model_name,
  updateCounts,
  counts,
  template,
  setAttribute,
  isAdminUser
}) => {
  const [showHiddenAttributes, setShowHiddenAttributes] = useState(false);
  const [canReparent, setCanReparent] = useState(false);
  const [determiningCanReparent, setDeterminingCanReparent] = useState(true);
  const dispatch = useDispatch();
  const invoke = useActionInvoker();
  const {executeAction} = useMagmaActions();

  const models = useReduxState((state) => selectModels(state));
  const classes = reportStyles();

  const modelCount = counts[model_name]?.count;
  const attributeCounts = counts[model_name]?.attributes;

  const getAnswer = (query, handle) =>
    requestAnswer({query})(dispatch).then(({answer}) => handle(answer));

  const countAttributes = () => {
    if (attributeCounts != undefined) return;

    Object.keys(template.attributes).forEach((attribute_name) => {
      let query = [model_name, ['::has', attribute_name], '::count'];

      getAnswer(query, (count) =>
        updateCounts({
          type: 'ATTRIBUTE_COUNT',
          model_name,
          attribute_name,
          count
        })
      );
    });
  };

  const [showDiff, setShowDiff] = useState(false);
  const [diffProject, setDiffProject] = useState(null);
  const [diffModel, setDiffModel] = useState(null);
  const [diffTemplate, setDiffTemplate] = useState(null);

  const [anchor, setAnchor] = useState(null);

  const getDiffModelTemplate = (project_name, model_name) => {
    setDiffProject(project_name);
    setDiffModel(model_name);
    getDocuments(
      {
        project_name,
        model_name,
        record_names: [],
        attribute_names: 'all'
      },
      fetch
    ).then(({models}) => setDiffTemplate(models[model_name].template));
  };

  const handleAddAttribute = useCallback(
    (params) => {
      executeAction(
        addAttributes([{
          model_name,
          ...params
        }])
      );
    },
    [model_name, executeAction]
  );

  const handleAddLink = useCallback(
    (params) => {
      executeAction(
        addLink({
          modelName: model_name,
          ...params
        })
      );
    },
    [model_name]
  );

  const handleRemoveModel = useCallback(() => {
    executeAction(removeModel({model_name})).then(() =>
      invoke(removeModelAction(model_name))
    );
  }, [model_name]);

  const handleReparentModel = useCallback(
    (parent_model_name) => {
      executeAction(reparentModel({model_name, parent_model_name}));
    },
    [model_name]
  );

  const handleCopyModel = useCallback(
    (attributes) => {
      executeAction(
        addAttributes(attributes)
      );
    },
    [model_name]
  );

  const handleAddModel = useCallback(
    (params) => {
      executeAction(
        addModel({
          ...params,
          parent_model_name: model_name
        })
      );
    },
    [model_name]
  );

  const isLeaf = useMemo(() => {
    if (!models || !model_name || !models[model_name]) return;

    return isLeafModel(models[model_name]);
  }, [model_name, models]);

  useEffect(() => {
    if (counts[model_name]?.count >= 0) {
      setDeterminingCanReparent(false)
      setCanReparent(counts[model_name].count <= 0);
    } else {
      getAnswer([model_name, '::count'], (count) => {
        updateCounts({type: 'MODEL_COUNT', model_name, count});
        setDeterminingCanReparent(false)
        setCanReparent(0 == count);
      });
    }
  }, [model_name, counts]);

  return (
    <Grid className={classes.model_report}>
      <MapHeading className={classes.heading} name='Model' title={model_name}>
        {diffTemplate && (
          <Chip
            label={`diff: ${diffProject}.${diffModel}`}
            onDelete={() => {
              setDiffTemplate(null);
              setDiffModel(null);
              setDiffProject(null);
            }}
          />
        )}
        {modelCount != undefined && modelCount >= 0 && (
          <Chip
            label={`${modelCount} ${
              modelCount > 1 || modelCount == 0 ? 'records' : 'record'
            }`}
          />
        )}
        <IconButton size='small' onClick={(e) => setAnchor(e.target)}>
          <MoreVertIcon />
        </IconButton>
        <Menu
          elevation={1}
          style={{marginTop: '40px'}}
          anchorEl={anchor}
          open={Boolean(anchor)}
          onClose={() => setAnchor(null)}
        >
          <MenuItem
            onClick={() => {
              setShowDiff(true);
              setAnchor(null);
            }}
          >
            Compare with another model
          </MenuItem>
          <MenuItem
            disabled={attributeCounts != undefined}
            onClick={() => {
              countAttributes(model_name);
              setAnchor(null);
            }}
          >
            Count attributes
          </MenuItem>
          {isAdminUser && (
            <MenuItem
              onClick={() => {
                setShowHiddenAttributes(!showHiddenAttributes);
                setAnchor(null);
              }}
            >
              Toggle hidden attributes
            </MenuItem>
          )}
        </Menu>
      </MapHeading>
      {isAdminUser && (
        <ManageModelActions
          handleAddAttribute={handleAddAttribute}
          handleAddLink={handleAddLink}
          handleRemoveModel={handleRemoveModel}
          handleReparentModel={handleReparentModel}
          handleAddModel={handleAddModel}
          handleCopyModel={handleCopyModel}
          isLeaf={isLeaf}
          canReparent={canReparent}
          determiningCanReparent={determiningCanReparent}
          modelName={model_name}
        />
      )}
      <SelectProjectModelDialog
        open={showDiff}
        onClose={() => setShowDiff(false)}
        update={(project_name, model_name) =>
          getDiffModelTemplate(project_name, model_name)
        }
        title='Compare Models'
        buttonLabel='Compare'
        description='Select a comparison project and model'
      />
      <ModelAttributesTable
        className={classes.attributes}
        template={template}
        diffTemplate={diffTemplate}
        attributeCounts={attributeCounts}
        showHiddenAttributes={showHiddenAttributes}
        setAttribute={setAttribute}/>
    </Grid>
  );
};

export default ModelReport;
