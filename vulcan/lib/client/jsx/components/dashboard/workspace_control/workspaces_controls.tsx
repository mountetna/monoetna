import React, {
  useCallback,
  useContext,
  useState,
  useEffect,
  useMemo
} from 'react';
import 'regenerator-runtime/runtime';

import Grid from '@material-ui/core/Grid';
import ImageList from '@material-ui/core/ImageList';
import ImageListItem from '@material-ui/core/ImageListItem';
import TextField from '@material-ui/core/TextField';
import {makeStyles} from '@material-ui/core/styles';
import InputAdornment from '@material-ui/core/InputAdornment';
import SearchIcon from '@material-ui/icons/Search';
import Autocomplete from '@material-ui/lab/Autocomplete';

import {VulcanContext} from '../../../contexts/vulcan_context';
import {Workspaces} from '../../../api_types';
import Tag from '../tag';

const useStyles = makeStyles((theme) => ({
  controls: {
    padding: '15px',
    width: '100%',
    flex: '1 1 auto'
  },
  tagauto: {},
  tags: {
    padding: '12.5px !important'
  }
}));

export default function WorkspacesControls({
  project_name,
  setSearchString,
  setTags,
  tags,
  searchString
}: {
  project_name: string;
  setSearchString: Function;
  setTags: Function;
  tags?: string[];
  searchString?: string;
}) {
  const {showErrors, getWorkspaces} = useContext(VulcanContext);

  const [allWorkspaces, setAllWorkspaces] = useState<Workspaces>([]);

  useEffect(() => {
    showErrors(
      getWorkspaces(
        project_name
      ).then(({workspaces}: {workspaces: Workspaces}) =>
        setAllWorkspaces(workspaces)
      )
    );
  }, []);

  const classes = useStyles();

  const allTags = useMemo(() => {
    return [
      ...new Set(
        allWorkspaces.reduce(
          (acc: string[], f) => acc.concat(f.tags || []),
          ['public']
        )
      )
    ];
  }, [allWorkspaces]);

  return (
    <Grid item container className={classes.controls} spacing={6}>
      <Grid item xs={6}>
        <TextField
          fullWidth
          label='Search'
          variant='outlined'
          InputProps={{
            startAdornment: (
              <InputAdornment position='start'>
                <SearchIcon />
              </InputAdornment>
            )
          }}
          value={searchString}
          onChange={(e) => setSearchString(e.target.value)}
        />
      </Grid>
      <Grid item xs={6}>
        <Autocomplete
          fullWidth
          multiple
          className='figure-tag-autocomplete'
          classes={{
            input: classes.tags
          }}
          value={tags}
          options={allTags.sort()}
          renderInput={(params: any) => (
            <TextField
              {...params}
              size='small'
              label='Tags'
              variant='outlined'
            />
          )}
          renderOption={(option, state) => <span>{option}</span>}
          renderTags={(tags, getTagProps) =>
            tags.map((tag, index) => (
              <Tag {...getTagProps({index})} label={tag} />
            ))
          }
          filterOptions={(options, state) => {
            let regex = new RegExp(state.inputValue);
            return options.filter((o) => regex.test(o));
          }}
          onChange={(e, v) => setTags(v)}
        />
      </Grid>
    </Grid>
  );
}
