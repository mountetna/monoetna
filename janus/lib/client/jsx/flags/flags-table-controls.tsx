import React, {useState, useEffect, useCallback} from 'react';
import Grid from '@material-ui/core/Grid';
import Button from '@material-ui/core/Button';
import {makeStyles} from '@material-ui/core/styles';
import TextField from '@material-ui/core/TextField';
import InputAdornment from '@material-ui/core/InputAdornment';
import Search from '@material-ui/icons/Search';
import Input from '@material-ui/core/Input';
import InputLabel from '@material-ui/core/InputLabel';
import MenuItem from '@material-ui/core/MenuItem';
import FormControl from '@material-ui/core/FormControl';
import Select from '@material-ui/core/Select';
import Chip from '@material-ui/core/Chip';
import AddIcon from '@material-ui/icons/Add';
import RemoveIcon from '@material-ui/icons/Remove';
import IconButton from '@material-ui/core/IconButton';
import Tooltip from '@material-ui/core/Tooltip';
import Autocomplete from '@material-ui/lab/Autocomplete';

import {fetchProjects} from '../api/janus_api';
import {Project} from '../types/janus_types';

const useStyles = makeStyles((theme) => ({
  chips: {
    display: 'flex',
    flexWrap: 'wrap'
  },
  formControl: {
    margin: theme.spacing(1),
    minWidth: 300,
    maxWidth: 400
  },
  chip: {
    margin: 2
  },
  button: {
    margin: theme.spacing(1)
  }
}));

// Assumes this won't be in any project or flag string...
const DELIMITER = ' -- ';

function MultiSelector({
  options,
  label,
  onChange
}: {
  options: string[];
  label: string;
  onChange: (selection: string[]) => void;
}) {
  const classes = useStyles();

  function onSelect(selection: string[]) {
    onChange(selection.map((s) => unpackOption(s)));
  }

  function unpackOption(opt: string) {
    return opt.split(DELIMITER)[0].trim();
  }

  return (
    <FormControl className={classes.formControl}>
      <Autocomplete
        multiple
        id={`${label}-filter`}
        options={options}
        getOptionLabel={(option) => unpackOption(option)} // This shows up in the chip
        onChange={(e, value) => onSelect(value as string[])}
        renderInput={(params: any) => (
          <TextField
            {...params}
            variant='outlined'
            label={label}
            placeholder={label}
          />
        )}
        renderOption={(option, state) => <span>{option}</span>}
        filterOptions={(options, state) => {
          let regex = new RegExp(state.inputValue);
          return options.filter((o) => regex.test(o));
        }}
      />
    </FormControl>
  );
}

const TableControls = ({
  onChangeSearch,
  onChangeProjects,
  onChangeFlags,
  flagOptions,
  projectOptions
}: {
  flagOptions: string[];
  projectOptions: string[];
  onChangeSearch: (search: string) => void;
  onChangeProjects: (projects: string[]) => void;
  onChangeFlags: (flags: string[]) => void;
}) => {
  const [projects, setProjects] = useState([] as Project[]);
  const [prettyProjectOptions, setPrettyProjectOptions] = useState(
    [] as string[]
  );
  const classes = useStyles();

  useEffect(() => {
    fetchProjects().then(({projects}) => setProjects(projects));
  }, []);

  useEffect(() => {
    setPrettyProjectOptions(
      projectOptions.map((p: string): string => {
        let fullProject = projects.find(
          (project) => project.project_name === p
        );

        return fullProject
          ? `${p}${DELIMITER}${fullProject.project_name_full}`
          : p;
      })
    );
  }, [projectOptions, projects]);

  return (
    <Grid container xs={12}>
      <Grid item xs={3}>
        <TextField
          label='Search'
          variant='outlined'
          onChange={(e) => onChangeSearch(e.target.value as string)}
          InputLabelProps={{
            shrink: true
          }}
          InputProps={{
            startAdornment: (
              <InputAdornment position='start'>
                <Search />
              </InputAdornment>
            )
          }}
        />
      </Grid>
      <Grid item xs={3}>
        <MultiSelector
          options={prettyProjectOptions}
          onChange={onChangeProjects}
          label='Projects'
        />
      </Grid>
      <Grid item xs={3}>
        <MultiSelector
          options={flagOptions}
          onChange={onChangeFlags}
          label='Flags'
        />
      </Grid>
      <Grid item xs={3}>
        <Button
          variant='contained'
          color='primary'
          className={classes.button}
          startIcon={<AddIcon />}
        >
          Add Flag
        </Button>
        <Button
          variant='contained'
          color='secondary'
          className={classes.button}
          startIcon={<RemoveIcon />}
        >
          Remove Flag
        </Button>
      </Grid>
    </Grid>
  );
};

export default TableControls;
