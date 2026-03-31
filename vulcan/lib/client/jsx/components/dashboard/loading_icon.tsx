import React from 'react';
import AutorenewIcon from '@material-ui/icons/Autorenew';
import {makeStyles} from '@material-ui/core/styles';
import Typography from '@material-ui/core/Typography';
import Tooltip from '@material-ui/core/Tooltip';
import Grid from '@material-ui/core/Grid';
import Icon from 'etna-js/components/icon';

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

export function LoadingIconWithText({text=''}: {text:string}) {
    return <Grid container>
        <Grid item>
            <LoadingIcon/>
        </Grid>
        <Grid item>
            <Typography>{text}</Typography>
        </Grid>
    </Grid>
};

export function FlatLoadbleIcon({label, tooltip, loading=false, disabled=false}: {
    label: string;
    tooltip: string;
    loading: boolean;
    disabled; boolean;
}) {
    const classes = useStyles();
    return <div title={tooltip}>
        <Grid container>
            <Grid item>
                { label && <span className='flat-label'>{label}</span> }
            </Grid>
            <Grid item>
                    { loading ? <Icon disabled={disabled} icon={'spinner fa-spin'} /> : <AutorenewIcon color={'disabled'}/> }
            </Grid>
        </Grid>
    </div>
};
