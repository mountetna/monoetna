import React, { useState, useEffect, useCallback, useContext } from 'react';
import Typography from '@material-ui/core/Typography';
import { makeStyles } from '@material-ui/core/styles';
import Grid from '@material-ui/core/Grid';
import Button from '@material-ui/core/Button';
import TextField from '@material-ui/core/TextField';
import FormControl from '@material-ui/core/FormControl';
import Select from '@material-ui/core/Select';
import MenuItem from '@material-ui/core/MenuItem';

import ProjectHeader from 'etna-js/components/project-header';
import { json_get } from 'etna-js/utils/fetch';
import { magmaPath } from 'etna-js/api/magma_api';

import { MagmaRulesResponse } from './utils/rules';



const useStyles = makeStyles((theme) => ({
  frame: {
    height: 'calc(100vh - 61px - 48px)'
  },
  main: {
    width: 'auto',
    flex: '1 1 auto',
    padding: '100px 0px 200px',
  },
  admin: {
    flex: '1 1 auto',
    width: 'auto',
    borderLeft: '1px solid #eee',
    margin: '20px 0px'
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

const GnomonMain = ({ project_name }: { project_name: string }) => {
  const classes = useStyles();

  const [identifier, setIdentifier] = useState('');

  const [rules, setRules] = useState<string[] | null>(null);

  const onEnterIdentifier = useCallback(e => {
    if (e.key == 'Enter') {
      window.location.href = `/${project_name}/identify/${identifier}`;
    }
  }, [identifier]);

  const onSelectRule = useCallback(e => {
    window.location.href = `/${project_name}/create/${e.target.value}`;
  }, []);

  const onClickBulkCreate = () => {
    window.location.href = `/${project_name}/create`
  }

  const onClickBrowse = () => {
    window.location.href = `/${project_name}/browse`
  }

  useEffect(() => {
    async function fetchRules() {
      try {
        const { config }: { config: MagmaRulesResponse } = await json_get(magmaPath(`gnomon/${project_name}`))
        setRules(Object.keys(config.rules))
      } catch (error) {
        setRules([])
        console.error(`Error fetching rules: ${error}`)
      }
    }

    fetchRules()
  }, []);

  return (
    <Grid>
      <ProjectHeader project_name={project_name} className={classes.header} />

      <Grid container className={classes.frame} direction='row'>

        <Grid item className={classes.main} alignItems='center' container justifyContent='space-around' direction='column'>
          <Grid>
            <TextField variant='outlined'
              onKeyPress={onEnterIdentifier}
              className={classes.text}
              placeholder='Enter an identifier'
              value={identifier}
              onChange={(e) => setIdentifier(e.target.value)}
            />
          </Grid>

          <Grid><Typography color='secondary'>∼ OR ∼</Typography></Grid>

          <Grid>
            <Button
              color="secondary"
              size="large"
              onClick={onClickBrowse}
            >
              Browse identifiers
            </Button>
          </Grid>

          <Grid><Typography color='secondary'>∼ OR ∼</Typography></Grid>

          <Grid>
            <FormControl variant='outlined'>
              <Select value='' onChange={onSelectRule} displayEmpty>
                <MenuItem value='' disabled>Create an identifier</MenuItem>
                {
                  rules && rules?.map(rule => <MenuItem key={rule} value={rule}>{rule}</MenuItem>)
                }
              </Select>
            </FormControl>
          </Grid>

          <Grid><Typography color='secondary'>∼ OR ∼</Typography></Grid>

          <Grid>
            <Button
              color="secondary"
              size="large"
              onClick={onClickBulkCreate}
            >
              Create many identifiers
            </Button>
          </Grid>
        </Grid>

        <Grid item container alignItems='center' justifyContent='space-around' className={classes.admin}>
          <Button color='secondary' size='large'
            onClick={() => window.location.href = `/${project_name}/rules/`}>Edit Rules</Button>
        </Grid>

      </Grid>

    </Grid>
  );
};

export default GnomonMain;
