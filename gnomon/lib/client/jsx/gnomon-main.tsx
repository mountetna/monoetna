import React, {useState, useEffect, useCallback, useContext} from 'react';
import { useDispatch } from 'react-redux';
import Typography from '@material-ui/core/Typography';
import {makeStyles} from '@material-ui/core/styles';
import Grid from '@material-ui/core/Grid';
import Button from '@material-ui/core/Button';
import TextField from '@material-ui/core/TextField';
import FormControl from '@material-ui/core/FormControl';
import Select from '@material-ui/core/Select';
import MenuItem from '@material-ui/core/MenuItem';
import ProjectHeader from 'etna-js/components/project-header';

import {requestAnswer} from 'etna-js/actions/magma_actions';

const useStyles = makeStyles((theme) => ({
  main: {
    width: '100%',
    padding: '100px 0px 200px',
    height: 'calc(100vh - 61px)'
  },
  header: {
    borderBottom: '1px solid #eee'
  },
  text: {
    fontFamily: 'monospaced',
    '& input::placeholder': {
      fontFamily: 'Open Sans, sans-serif',
      textAlign: 'center'
    }
  }
}));

const GnomonMain = ({project_name}: {project_name: string}) => {
  const classes = useStyles();

  const [ identifier, setIdentifier ] = useState('');

  const [ models, setModels ] = useState(null);
  const [ error, setError ] = useState(null);

  const onEnter = useCallback( e => {
    if (e.key == 'Enter') {
      window.location.href = `/${project_name}/identify/${identifier}`;
    }
  }, [ identifier ]);

  const dispatch = useDispatch();

  useEffect(() => {
    requestAnswer({
      project_name,
      query: '::model_names'
    })(dispatch).then(
      ({answer}) => { setModels(answer.sort()); setError(null); }
    ).catch( e => e.then( ({error}) => setError(error) ) );
  }, []);

  return (
    <Grid>
      <ProjectHeader project_name={project_name} className={classes.header}/>
      <Grid className={classes.main} alignItems='center' container justify='space-around' direction='column'>
        <Grid>
          <TextField variant='outlined'
            onKeyPress={ onEnter }
            className={classes.text}
            placeholder='Enter an identifier'
            value={ identifier }
            onChange={ (e) => setIdentifier(e.target.value) }
            />
        </Grid>
        <Grid><Typography color='secondary'>∼ OR ∼</Typography></Grid>
        <Grid>
          <FormControl variant='outlined'>
            <Select value='' displayEmpty>
              <MenuItem value='' disabled>Create an identifier</MenuItem>
              {
                models && models.map( m => <MenuItem key={m} value={m}>{m}</MenuItem> )
              }
            </Select>
          </FormControl>
        </Grid>
      </Grid>
    </Grid>
  );
};

export default GnomonMain;
