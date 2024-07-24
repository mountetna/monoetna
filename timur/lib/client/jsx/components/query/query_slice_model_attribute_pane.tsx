import React, {useContext, useState, useMemo} from 'react';
import Grid from '@material-ui/core/Grid';
import IconButton from '@material-ui/core/IconButton';
import AddIcon from '@material-ui/icons/Add';
import Tooltip from '@material-ui/core/Tooltip';
import Typography from '@material-ui/core/Typography';

import {makeStyles} from '@material-ui/core/styles';

import {QueryGraphContext} from '../../contexts/query/query_graph_context';
import QuerySliceControl from './query_slice_control';
import {QuerySlice, QueryColumn} from '../../contexts/query/query_types';
import useSliceMethods from './query_use_slice_methods';

const useStyles = makeStyles((theme) => ({
  slices: {
    paddingLeft: '25px'
  }
}));

const QuerySliceModelAttributePane = ({
  column,
  columnIndex,
  showControls
}: {
  column: QueryColumn;
  columnIndex: number;
  showControls: boolean;
}) => {
  // All the slices related to a given model / attribute,
  //   with the model / attribute as a "label".
  // Matrices will have modelName + attributeName.
  const [updateCounter, setUpdateCounter] = useState(0);
  const {
    state: {graph, rootModel}
  } = useContext(QueryGraphContext);
  const classes = useStyles();

  const {matrixModelNames, addNewSlice, handlePatchSlice, handleRemoveSlice} =
    useSliceMethods(columnIndex, updateCounter, setUpdateCounter);

  let sliceableModelNames = useMemo(() => {
    if (!rootModel) return [];

    if (
      rootModel === column.model_name &&
      matrixModelNames.includes(column.model_name)
    ) {
      return [column.model_name];
    } else {
      return graph
        .sliceableModelNamesInPath(rootModel, column.model_name)
        .sort();
    }
  }, [column, graph, rootModel, matrixModelNames]);

  return (
    <Grid container className={classes.slices}>
      <Grid item container alignItems='center'>
        Where:
        { true && <Tooltip title='Add slice' aria-label='Add slice'>
            <IconButton size='small' onClick={addNewSlice} color='primary'>
              <AddIcon fontSize='small'/>
            </IconButton>
          </Tooltip>
        }
      </Grid>
      {
      (column && column.slices.length) ?
      <Grid container direction='column'>
        {column?.slices.map((slice: QuerySlice, index: number) => (
            <QuerySliceControl
              key={`model-${column.model_name}-${index}-${updateCounter}`}
              slice={slice}
              sliceIndex={index}
              modelNames={sliceableModelNames}
              graph={graph}
              patchSlice={(updatedSlice: QuerySlice) =>
                handlePatchSlice(index, updatedSlice)
              }
              removeSlice={() => handleRemoveSlice(index)}
            />
        ))}
<<<<<<< HEAD
      </Grid>
      <Tooltip title='Add slice' aria-label='Add slice'>
        <Button
          aria-label='Add slice'
          onClick={() => addNewSlice()}
          startIcon={<AddIcon />}
        >
          Add slice
        </Button>
      </Tooltip>
    </React.Fragment>
=======
      </Grid> : <Typography component='div' color='gray'>no conditions</Typography>
      }
    </Grid>
>>>>>>> ec64dc12d (fix select columns)
  );
};

export default QuerySliceModelAttributePane;
