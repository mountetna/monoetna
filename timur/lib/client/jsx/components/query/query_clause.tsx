import React, {PropsWithChildren} from 'react';
import Grid from '@material-ui/core/Grid';
import {makeStyles} from '@material-ui/core/styles';

const useStyles = makeStyles((theme) => ({
  clauseTitle: {
    fontSize: '1.2rem'
  },
  queryClause: {
    padding: '15px',
    paddingLeft: '5px',
    borderBottom: '1px solid #ccc'
  }
}));

const QueryClause = ({title, children}: PropsWithChildren<{title: string}>) => {
  const classes = useStyles();

  return (
    <Grid alignItems='flex-start' container className={classes.queryClause}>
      <Grid xs={1} item className={classes.clauseTitle}>
        {title}
      </Grid>
      <Grid xs={11} item>
        {children}
      </Grid>
    </Grid>
  );
};

export default QueryClause;
