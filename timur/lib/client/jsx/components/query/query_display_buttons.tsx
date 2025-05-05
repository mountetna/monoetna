import React from 'react';

import Grid from '@material-ui/core/Grid';
import Tooltip from '@material-ui/core/Tooltip';
import ToggleButton from '@material-ui/lab/ToggleButton';
import ToggleButtonGroup from '@material-ui/lab/ToggleButtonGroup';
import ListIcon from '@material-ui/icons/List';
import CodeIcon from '@material-ui/icons/Code';
import TuneIcon from '@material-ui/icons/Tune';
import SearchHistoryIcon from '@material-ui/icons/YoutubeSearchedFor';

import { makeStyles } from '@material-ui/core/styles';

import QueryOptions from './query_options';
import QueryHistory from './query_history';

const useStyles = makeStyles((theme) => ({
  buttons: {
    borderRadius: '2px'
  }
}));


const QueryDisplayButtons = ({options,setOptions}:{
  options: string[];
  setOptions: (opts: string[]) => void;
}) => {
  const classes = useStyles();

  const showQuery = options.includes('controls');
  const showQueryOptions = options.includes('options');
  const showQueryHistory = options.includes('history');

  return (
    <>
      <ToggleButtonGroup className={classes.buttons} value={options} onChange={ (e,f) => setOptions(f) } size='small'>
        <ToggleButton aria-label="controls" disableRipple={true} value='controls'>
          <Tooltip title="Query controls">
            <ListIcon/>
          </Tooltip>
        </ToggleButton>
        { showQuery &&
          <ToggleButton aria-label="raw" disableRipple={true} value='raw'>
            <Tooltip title="Raw query">
              <CodeIcon/>
            </Tooltip>
          </ToggleButton>
        }
        <ToggleButton aria-label="options" disableRipple={true} value='options'>
          <Tooltip title="Query Options">
            <TuneIcon/>
          </Tooltip>
        </ToggleButton>
        <ToggleButton aria-label="history" disableRipple={true} value='history'>
          <Tooltip title="Search History">
            <SearchHistoryIcon/>
          </Tooltip>
        </ToggleButton>
      </ToggleButtonGroup>
      <QueryOptions open={showQueryOptions} onClose={ () => setOptions(options.filter( o => o !== 'options')) }/>
      <QueryHistory open={showQueryHistory} onClose={ () => setOptions(options.filter( o => o !== 'history')) }/>
    </>
  );
};

export default QueryDisplayButtons;
