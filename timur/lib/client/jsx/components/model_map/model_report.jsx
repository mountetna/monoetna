import {useDispatch} from 'react-redux';
import React, {useState, useCallback, useMemo, useEffect} from 'react';
import {sortAttributeList} from '../../utils/attributes';
import SelectProjectModelDialog from '../select_project_model';
import {requestAnswer} from 'etna-js/actions/magma_actions';
import {getDocuments} from 'etna-js/api/magma_api';
import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {isEqual} from 'lodash';

import Grid from '@material-ui/core/Grid';
import Typography from '@material-ui/core/Typography';
import TextField from '@material-ui/core/TextField';
import {makeStyles} from '@material-ui/core/styles';
import Menu from '@material-ui/core/Menu';
import MenuItem from '@material-ui/core/MenuItem';
import Table from '@material-ui/core/Table';
import TableBody from '@material-ui/core/TableBody';
import TableCell from '@material-ui/core/TableCell';
import TableSortLabel from '@material-ui/core/TableSortLabel';
import TableContainer from '@material-ui/core/TableContainer';
import TableHead from '@material-ui/core/TableHead';
import TableRow from '@material-ui/core/TableRow';
import Paper from '@material-ui/core/Paper';
import MapHeading from './map_heading';
import Button from '@material-ui/core/Button';
import InputAdornment from '@material-ui/core/InputAdornment';
import IconButton from '@material-ui/core/IconButton';
import SearchIcon from '@material-ui/icons/Search';
import MoreVertIcon from '@material-ui/icons/MoreVert';
import Tooltip from '@material-ui/core/Tooltip';
import LinearProgress from '@material-ui/core/LinearProgress';
import Chip from '@material-ui/core/Chip';
import AddIcon from '@material-ui/icons/Add';
import LinkIcon from '@material-ui/icons/Link';
import DeleteIcon from '@material-ui/icons/Delete';
import SwapHorizIcon from '@material-ui/icons/SwapHoriz';
import LibraryAddIcon from '@material-ui/icons/LibraryAdd';
import VisibilityOffIcon from '@material-ui/icons/VisibilityOff';

import AddAttributeModal from './add_attribute_modal';
import AddLinkModal from './add_link_modal';
import ReparentModelModal from './reparent_model_modal';
import AddModelModal from './add_model_modal';
import RemoveModelModal from './remove_model_modal';
import {useModal} from 'etna-js/components/ModalDialogContainer';
import {selectModels} from 'etna-js/selectors/magma';
import {useReduxState} from 'etna-js/hooks/useReduxState';
import {
  addAttribute,
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

const attributeStyles = makeStyles((theme) => ({
  attribute: {
    wordBreak: 'break-all'
  },
  value: {
    color: 'darkgoldenrod',
    cursor: 'pointer'
  },
  missing: {},
  indicator: {
    width: '30px',
    color: 'gray',
    cursor: 'default'
  },
  addBtn: {
    margin: '0.5rem'
  },
  add: {
    width: '30px',
    color: 'green'
  },
  type: {
    color: 'gray',
    width: '25%',
    paddingRight: '10px'
  },
  progress: {
    width: '30px',
    margin: '5px'
  },
  counts: {
    width: '10%',
    minWidth: '120px'
  },
  ident: {},
  changed: {
    backgroundColor: 'rgba(255,255,0,0.1)'
  },
  present: {
    backgroundColor: 'rgba(0,255,0,0.1)'
  },
  absent: {
    backgroundColor: 'rgba(255,0,0,0.1)'
  },
  hiddenIcon: {
    marginRight: '1rem'
  },
  typeWrapper: {
    display: 'flex',
    justifyContent: 'right'
  },
  hiddenTypeWrapper: {
    display: 'flex',
    justifyContent: 'right',
    paddingTop: '0.35rem'
  }
}));

const diffTypes = {
  ident: {ind: '', title: ''},
  present: {ind: '+', title: 'Present in this model'},
  absent: {ind: '-', title: 'Absent in this model'},
  changed: {ind: 'c', title: 'Changed in this model'}
};

const ManageModelActions = ({
  handleAddAttribute,
  handleAddLink,
  handleAddModel,
  handleRemoveModel,
  handleReparentModel,
  isLeaf,
  canReparent,
  modelName
}) => {
  const classes = attributeStyles();
  const {openModal} = useModal();

  return (
    <>
      <Tooltip title='Add Attribute' aria-label='Add Attribute'>
        <Button
          className={classes.addBtn}
          startIcon={<AddIcon />}
          onClick={() => {
            openModal(<AddAttributeModal onSave={handleAddAttribute} />, {
              closeOnClickBackdrop: false
            });
          }}
        >
          Attribute
        </Button>
      </Tooltip>
      <Tooltip title='Add Link' aria-label='Add Link'>
        <Button
          className={classes.addBtn}
          startIcon={<LinkIcon />}
          onClick={() => {
            openModal(<AddLinkModal onSave={handleAddLink} />, {
              closeOnClickBackdrop: false
            });
          }}
        >
          Link
        </Button>
      </Tooltip>
      <Tooltip title='Add Model' aria-label='Add Model'>
        <Button
          className={classes.addBtn}
          startIcon={<LibraryAddIcon />}
          onClick={() => {
            openModal(
              <AddModelModal modelName={modelName} onSave={handleAddModel} />,
              {
                closeOnClickBackdrop: false
              }
            );
          }}
        >
          Model
        </Button>
      </Tooltip>
      <Tooltip title='Reparent Model' aria-label='Reparent Model'>
        <Button
          className={classes.addBtn}
          startIcon={<SwapHorizIcon />}
          disabled={!canReparent}
          onClick={() => {
            openModal(
              <ReparentModelModal
                modelName={modelName}
                onSave={handleReparentModel}
              />,
              {
                closeOnClickBackdrop: false
              }
            );
          }}
        >
          Reparent Model
        </Button>
      </Tooltip>
      {isLeaf && (
        <Tooltip title='Remove Model' aria-label='Remove Model'>
          <Button
            className={classes.addBtn}
            startIcon={<DeleteIcon />}
            onClick={() => {
              openModal(
                <RemoveModelModal
                  modelName={modelName}
                  onSave={handleRemoveModel}
                />,
                {
                  closeOnClickBackdrop: false
                }
              );
            }}
          >
            Remove Model
          </Button>
        </Tooltip>
      )}
    </>
  );
};

const ModelAttribute = ({
  attribute_name,
  template,
  diffTemplate,
  setAttribute,
  count,
  modelCount,
  isHidden
}) => {
  const classes = attributeStyles();

  const attribute = template?.attributes[attribute_name];

  const diffAttribute = diffTemplate?.attributes[attribute_name];

  const [displayAttribute, diffType] = !diffTemplate
    ? [attribute, 'ident']
    : attribute && diffAttribute
    ? isEqual(attribute, diffAttribute)
      ? [attribute, 'ident']
      : [attribute, 'changed']
    : attribute
    ? [attribute, 'present']
    : [diffAttribute, 'absent'];

  const {attribute_type, attribute_group, description} = displayAttribute;

  if (!template) return null;

  return (
    <TableRow className={`${classes.attribute} ${classes[diffType]}`}>
      {diffTemplate && (
        <TableCell
          className={classes.indicator}
          align='left'
          title={diffTypes[diffType].title}
        >
          {diffTypes[diffType].ind}
        </TableCell>
      )}
      <TableCell className={classes.type} align='right'>
        <div
          className={isHidden ? classes.hiddenTypeWrapper : classes.typeWrapper}
        >
          {isHidden ? (
            <div className={classes.hiddenIcon}>
              <VisibilityOffIcon />
            </div>
          ) : null}
          {attribute_type}
        </div>
      </TableCell>
      <TableCell
        className={attribute ? classes.value : classes.missing}
        align='left'
        onClick={attribute ? () => setAttribute(attribute_name) : undefined}
      >
        {attribute_name}{' '}
      </TableCell>
      <TableCell align='left'>{attribute_group}</TableCell>
      <TableCell align='left'>{description}</TableCell>
      {count != undefined && (
        <TableCell className={classes.counts} align='left'>
          <Grid container alignItems='center'>
            {count}
            <LinearProgress
              className={classes.progress}
              variant='determinate'
              value={(100 * count) / (modelCount || 1)}
            />
            <Typography variant='body2' color='textSecondary'>
              {Math.round((100 * count) / (modelCount || 1))}%
            </Typography>
          </Grid>
        </TableCell>
      )}
    </TableRow>
  );
};

const reportStyles = makeStyles((theme) => ({
  model_report: {
    flex: '1 1 50%',
    overflowY: 'scroll',
    padding: '15px 10px'
  },
  heading: {
    marginBottom: '10px'
  },
  report_row: {
    borderBottom: '1px solid #eee'
  },
  table: {
    boxSizing: 'border-box',
    height: 'calc(100% - 95px)',
    overflowY: 'scroll'
  },
  filter: {
    paddingBottom: '10px'
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

const ATT_KEYS = {
  attribute: 'attribute_name',
  type: 'attribute_type',
  group: 'attribute_group',
  description: 'description'
};

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

  const [filterString, setFilterString] = useState('');

  const filterMatch = new RegExp(
    `^(?:(${Object.keys(ATT_KEYS).join('|')}):)?(.*)$`
  );

  const matchesFilter = useCallback(
    (attribute) => {
      if (filterString == '') return true;

      let tokens = filterString
        .split(/\s/)
        .filter((_) => _)
        .map((token) => token.match(filterMatch).slice(1));

      return tokens.every(([column, token]) => {
        const tokenMatch = new RegExp(token, 'i');
        const values = column
          ? [attribute[ATT_KEYS[column]]]
          : Object.values(ATT_KEYS).map((k) => attribute[k]);

        return values.some((s) => s?.match(tokenMatch));
      });
    },
    [filterString]
  );

  const [order, setOrder] = React.useState('asc');
  const [orderBy, setOrderBy] = React.useState('type');

  const sortByOrder = (attributes) => {
    let srt;
    if (orderBy === 'type') srt = sortAttributeList(attributes);
    else {
      srt = attributes.sort((a, b) =>
        (a[ATT_KEYS[orderBy]] || '').localeCompare(b[ATT_KEYS[orderBy]] || '')
      );
    }
    return order === 'desc' ? srt.reverse() : srt;
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
  const attributes = Object.values({
    ...(template?.attributes || {}),
    ...(diffTemplate && diffTemplate.attributes)
  });

  const handleAddAttribute = useCallback(
    (params) => {
      executeAction(
        addAttribute({
          model_name,
          ...params
        })
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
      setCanReparent(counts[model_name].count <= 0);
    } else {
      getAnswer([model_name, '::count'], (count) => {
        updateCounts({type: 'MODEL_COUNT', model_name, count});
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
      <TextField
        fullWidth
        placeholder='Filter attributes, e.g. "rna type:file"'
        variant='outlined'
        size='small'
        className={classes.filter}
        value={filterString}
        InputProps={{
          startAdornment: (
            <InputAdornment position='start'>
              <SearchIcon />
            </InputAdornment>
          )
        }}
        onChange={(e) => setFilterString(e.target.value)}
      />
      <TableContainer component={Paper} className={classes.table}>
        <Table stickyHeader size='small'>
          <TableHead>
            <TableRow>
              {diffTemplate && <TableCell className={classes.indicator} />}
              {['type', 'attribute', 'group', 'description', 'counts'].map(
                (key) =>
                  (key !== 'counts' || attributeCounts) && (
                    <TableCell
                      key={key}
                      className={key in classes ? classes[key] : null}
                      align={key == 'type' ? 'right' : 'left'}
                    >
                      <TableSortLabel
                        active={orderBy === key}
                        direction={orderBy === key ? order : 'asc'}
                        onClick={((key) => (e) => {
                          setOrder(
                            orderBy === key && order === 'asc' ? 'desc' : 'asc'
                          );
                          setOrderBy(key);
                        })(key)}
                      >
                        {key}
                      </TableSortLabel>
                    </TableCell>
                  )
              )}
            </TableRow>
          </TableHead>

          <TableBody>
            {sortByOrder(attributes)
              .filter(
                (attribute) =>
                  (showHiddenAttributes ? true : !attribute.hidden) &&
                  matchesFilter(attribute)
              )
              .map((attribute) => (
                <ModelAttribute
                  key={attribute.attribute_name}
                  attribute_name={attribute.attribute_name}
                  setAttribute={setAttribute}
                  count={
                    attributeCounts && attributeCounts[attribute.attribute_name]
                  }
                  modelCount={modelCount}
                  template={template}
                  diffTemplate={diffTemplate}
                  isHidden={attribute.hidden}
                />
              ))}
          </TableBody>
        </Table>
        {isAdminUser && (
          <ManageModelActions
            handleAddAttribute={handleAddAttribute}
            handleAddLink={handleAddLink}
            handleRemoveModel={handleRemoveModel}
            handleReparentModel={handleReparentModel}
            handleAddModel={handleAddModel}
            isLeaf={isLeaf}
            canReparent={canReparent}
            modelName={model_name}
          />
        )}
      </TableContainer>
    </Grid>
  );
};

export default ModelReport;
