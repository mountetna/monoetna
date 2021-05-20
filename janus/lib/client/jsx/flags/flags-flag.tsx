import React from 'react';

import {makeStyles} from '@material-ui/core/styles';
import Typography from '@material-ui/core/Typography';
import LabelIcon from '@material-ui/icons/Label';

const useStyles = makeStyles({
  root: {
    position: 'relative',
    display: 'inline-flex',
    justifyContent: 'center',
    alignItems: 'center'
  },
  icon: {
    fontSize: '2.5em'
  },
  flag: {
    position: 'absolute',
    lineHeight: 1,
    color: '#fff',
    top: '0.5em',
    fontSize: '1em'
  }
});

const Flag = ({flag}: {flag: string}) => {
  const classes = useStyles();

  return (
    <div className={classes.root}>
      <LabelIcon color='primary' className={classes.icon} />
      <Typography component='span' className={classes.flag}>
        {flag}
      </Typography>
    </div>
  );
};

export default Flag;
