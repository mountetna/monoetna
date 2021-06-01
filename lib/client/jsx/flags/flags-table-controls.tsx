import React, {useState, useEffect, useContext} from 'react';
import Grid from '@material-ui/core/Grid';
import {makeStyles} from '@material-ui/core/styles';
import TextField from '@material-ui/core/TextField';
import InputAdornment from '@material-ui/core/InputAdornment';
import Search from '@material-ui/icons/Search';
import FormControl from '@material-ui/core/FormControl';
import Autocomplete from '@material-ui/lab/Autocomplete';

import {FlagsContext} from './flags-context';
import {fetchProjects} from '../api/janus_api';

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

  function onSelect(selection: string[]) {
    onChange(selection.map((s) => unpackOption(s)));
  }

  function unpackOption(opt: string) {
    return opt.split(DELIMITER)[0].trim();
  }

  return (
    <FormControl fullWidth={true}>
      <Autocomplete
        multiple
        id={`${label}-filter`}
        options={options.sort()}
        getOptionLabel={(option) => unpackOption(option)} // This shows up in the chip
        onChange={(e, value) => onSelect(value as string[])}
        renderInput={(params: any) => (
          <TextField
            {...params}
            label={label}
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
  const [prettyProjectOptions, setPrettyProjectOptions] = useState(
    [] as string[]
  );
  let {
    state: {projects},
    setProjects
  } = useContext(FlagsContext);

  useEffect(() => {
    fetchProjects().then(({projects}) => setProjects(projects));
  }, []);

  useEffect(() => {
    if (projects) {
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
    }
  }, [projectOptions, projects]);

  return (
    <React.Fragment>
      <Grid item xs={3}>
        <TextField
          label='Search'
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
    </React.Fragment>
  );
};

export default TableControls;
