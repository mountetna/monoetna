import React, {useState, useEffect, useCallback, useRef} from 'react';
import { json_get, json_post } from 'etna-js/utils/fetch';
import { copyText } from 'etna-js/utils/copy';
import {makeStyles} from '@material-ui/core/styles';

import moment from 'moment';
import MomentUtils from '@date-io/moment';
import TextField from '@material-ui/core/TextField';
import { DateTimePicker, MuiPickersUtilsProvider } from '@material-ui/pickers';
import EventIcon from '@material-ui/icons/Event';
import Button from '@material-ui/core/Button';
import Select from '@material-ui/core/Select';
import MenuItem from '@material-ui/core/MenuItem';
const useStyles = makeStyles((theme) => ({
  text: {
    width: '300px',
    marginRight: '5px'
  },
  date: {
    marginTop: '2px'
  }
}));


const TokenBuilder = ({user}) => {
  const [ payload, setPayload ] = useState(null);
  const [ date, setDate ] = useState(null);
  const [ error, setError ] = useState(null);

  let generateToken = useCallback(
    () => {
      let req;
      try {
        req = JSON.parse(payload);
      } catch (e) {
        setError(e.message)
        return;
      }

      req['exp'] = (new Date(date)).getTime();
      json_post('/api/tokens/build',req).then(
	(token) => { copyText(token); setError(null); }
      ).catch( p => p.then( ({error}) => setError(error) ) )
    }, [payload, date]
  );
  const classes = useStyles();

  return <div id='token-builder-group'>
    <div className='title'>Token Builder</div>
    { error && <div className='error'>{error}</div> }
    <div className='item'>
      <MuiPickersUtilsProvider libInstance={moment} utils={MomentUtils}>
        <TextField
          multiline
          className={classes.text}
          value={ payload }
          placeholder="JSON payload { email, perm }"
          onChange={ e => setPayload(e.target.value) }
          />
        <DateTimePicker className={classes.date} placeholder='expiration' leftArrowIcon={<EventIcon/>} clearable value={ date } onChange={ setDate }/>
        <Button
          onClick={ (payload && date) ? generateToken : null }
          disabled={ !payload || !date }>
          Copy Token
        </Button>
      </MuiPickersUtilsProvider>
    </div>
  </div>;
}

export default TokenBuilder;
