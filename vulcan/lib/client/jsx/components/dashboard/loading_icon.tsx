import React from 'react';
import AutorenewIcon from '@material-ui/icons/Autorenew';
import {makeStyles} from '@material-ui/core/styles';

const useStyles = makeStyles((theme) => ({
    loadingIcon: {
      'animation': '$spin 4s linear infinite'
    },
    '@keyframes spin': {
        '100%': {
            '-webkit-transform': 'rotate(360deg)',
            'transform': 'rotate(360deg)',
        }
    },
  }));

export default function LoadingIcon({}) {
    const classes = useStyles();
    return (
        <AutorenewIcon className={classes.loadingIcon}/>
    )
}
