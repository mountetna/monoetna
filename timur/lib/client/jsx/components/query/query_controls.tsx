import React, {useContext, useMemo} from 'react';
import Grid from '@material-ui/core/Grid';
import {makeStyles} from '@material-ui/core/styles';

import QueryRowPane from './query_row_pane';
import QueryColumnPane from './query_column_pane';
import QueryWherePane from './query_where_pane';

const useStyles = makeStyles((theme) => ({
  item: {
    width: '100%'
  }
}));

const QueryControls = () => {
  const classes = useStyles();

  return (
    <>
      <QueryRowPane />
      <QueryWherePane />
      <QueryColumnPane />
    </>
  );
};

export default QueryControls;
