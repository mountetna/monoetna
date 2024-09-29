import React, {useCallback, useState, useEffect, useMemo} from 'react';
import Grid from '@material-ui/core/Grid';
import FormControl from '@material-ui/core/FormControl';
import { makeStyles } from '@material-ui/core/styles';
import Autocomplete from '@material-ui/lab/Autocomplete';
import TextField from '@material-ui/core/TextField';
import Tooltip from '@material-ui/core/Tooltip';
import IconButton from '@material-ui/core/IconButton';
import FileCopyIcon from '@material-ui/icons/FileCopyOutlined';
import Typography from '@material-ui/core/Typography';
import Paper from '@material-ui/core/Paper';

import {Attribute} from '../../models/model_types';

import useSliceMethods from './query_use_slice_methods';
import {QueryColumn} from '../../contexts/query/query_types';
import {
  selectAllowedModelAttributes,
  attributeIsFile
} from '../../selectors/query_selector';
import QuerySliceModelAttributePane from './query_slice_model_attribute_pane';

import {visibleSortedAttributesWithUpdatedAt} from '../../utils/attributes';
import {QueryGraph} from '../../utils/query_graph';
import RemoveIcon from './query_remove_icon';
import Selector from './query_selector';
import QueryModelAttributeSelector from './query_model_attribute_selector';
import QueryNumber from './query_number';

const useStyles = makeStyles((theme) => ({
  fullWidth: {
    width: '80%',
    minWidth: 120
  },
  column_input: {
    flex: '1',
    '& input': {
      padding: '4px'
    }
  }
}));

function id(label: string) {
  return `${label}-${Math.random()}`;
}

const AttributeSelector = React.memo(
  ({
    onSelect,
    canEdit,
    attributeChoiceSet,
    column,
    label
  }: {
    label: string;
    column: QueryColumn;
    attributeChoiceSet: Attribute[];
    canEdit: boolean;
    onSelect: (attributeName: string) => void;
  }) => {
    const classes = useStyles();

    if (!canEdit)
      {return (
        <TextField
          variant="standard"
          disabled
          value={column.attribute_name}
          className='query-column-attribute' />
      );}

    return (
      <FormControl variant="standard" className={classes.fullWidth}>
        <Autocomplete
          id={`${id(label)}-attribute`}
          value={
            attributeChoiceSet.find(
              (a: Attribute) => a.attribute_name === column.attribute_name
            ) || null
          }
          options={attributeChoiceSet}
          getOptionLabel={(option) => option.attribute_name}
          renderInput={(params) => <TextField variant="standard" {...params} />}
          onChange={(e, v) => onSelect(v?.attribute_name || '')}
        />
      </FormControl>
    );
  }
);

const FilePredicateSelector = React.memo(
  ({
    onSelect,
    value
  }: {
    onSelect: (newValue: string) => void;
    value: string;
  }) => {
    return (
      <Selector
        canEdit={true}
        label='Predicate'
        name={value || 'url'}
        onSelect={onSelect}
        choiceSet={['url', 'md5', 'size']}
      />
    );
  }
);

const QueryColumnSelector = React.memo(
  ({
    label,
    column,
    columnIndex,
    canEdit,
    modelChoiceSet,
    graph,
    onSelectModel,
    onSelectAttribute,
    onChangeLabel,
    onRemoveColumn,
    onCopyColumn,
    onSelectPredicate
  }: {
    label: string;
    column: QueryColumn;
    columnIndex: number;
    canEdit: boolean;
    modelChoiceSet: string[];
    graph: QueryGraph;
    onSelectModel: (modelName: string) => void;
    onSelectAttribute: (attributeName: string) => void;
    onChangeLabel: (label: string) => void;
    onRemoveColumn: () => void;
    onCopyColumn: () => void;
    onSelectPredicate: (predicate: string) => void;
  }) => {
    const [selectableModelAttributes, setSelectableModelAttributes] = useState(
      [] as Attribute[]
    );
    // All the slices related to a given model / attribute,
    //   with the model / attribute as a "label".
    // Matrices will have modelName + attributeName.
    const [updateCounter, setUpdateCounter] = useState(0);

    const classes = useStyles();

    const selectAttributesForModel = useCallback(
      (modelName: string) => {
        let template = graph.template(modelName);
        setSelectableModelAttributes(
          selectAllowedModelAttributes(
            visibleSortedAttributesWithUpdatedAt(template.attributes)
          )
        );
      },
      [graph]
    );

    useEffect(() => {
      if (column?.model_name && graph.template(column.model_name)) {
        selectAttributesForModel(column.model_name);
      }
    }, [graph, column, selectAttributesForModel]);

    const {matrixModelNames, collectionModelNames} = useSliceMethods(
      columnIndex,
      updateCounter,
      setUpdateCounter
    );

    const isSliceableAsMatrix = matrixModelNames.includes(column.model_name);
    const isSliceableAsCollection = collectionModelNames.includes(
      column.model_name
    );

    const isSliceable = isSliceableAsMatrix || isSliceableAsCollection;

    const showFilePredicates = useMemo(() => {
      return (
        column?.attribute_name &&
        attributeIsFile(graph.models, column.model_name, column.attribute_name)
      );
    }, [column, graph]);

    const [ removeHint, setRemoveHint ] = useState(false);
    const [ showControls, setShowControls ] = useState(false);

    return (
      <Grid 
        style={{ textDecoration: removeHint ? 'line-through' : 'none' }}
        onMouseEnter={ () => setShowControls(true) }
        onMouseLeave={ () => setShowControls(false) }
        container direction='column'>
        <Grid
          item
          container
          alignItems='center'
          justifyContent='flex-start'
        >
          <QueryNumber
            setRemoveHint={ canEdit ? setRemoveHint : null }
            onClick={ canEdit ? onRemoveColumn : null }
            number={columnIndex} level={0}/>
          <QueryModelAttributeSelector
            setModel={ canEdit? onSelectModel : null }
            setAttribute={ canEdit ? onSelectAttribute : null}
            modelNames={modelChoiceSet}
            modelName={column.model_name}
            attributeName={column.attribute_name}
          />
          {
            column.model_name && 
              selectableModelAttributes.length > 0 &&
              showFilePredicates && <>
                &nbsp;
                <FilePredicateSelector
                  value={column.predicate || 'url'}
                  onSelect={onSelectPredicate}
                />
              </>
          }
          <Typography>&nbsp;as column&nbsp;</Typography>
            <TextField
              variant='standard'
              size='small'
              className={classes.column_input}
              value={column.display_label}
              onChange={(e) => onChangeLabel(e.target.value)}
            />
          { showControls && canEdit && <Tooltip title='Copy column' aria-label='Copy column'>
            <IconButton size='small' onClick={onCopyColumn} color='primary'>
              <FileCopyIcon fontSize='small'/>
            </IconButton>
          </Tooltip>}
        </Grid>
        {
          column.model_name && 
            selectableModelAttributes.length > 0 &&
            isSliceable &&
            canEdit &&
            <QuerySliceModelAttributePane showControls={ showControls } column={column} columnIndex={columnIndex} />
        }
      </Grid>
    );
  }
);

export default QueryColumnSelector;
