import React, {useState, useEffect, useMemo, useCallback} from 'react';
import {makeStyles} from '@material-ui/core/styles';
import DialogTitle from '@material-ui/core/DialogTitle';
import DialogContent from '@material-ui/core/DialogContent';
import DialogActions from '@material-ui/core/DialogActions';
import Grid from '@material-ui/core/Grid';
import TextField from '@material-ui/core/TextField';
import SearchIcon from '@material-ui/icons/Search';
import CancelIcon from '@material-ui/icons/Cancel';
import CircularProgress from '@material-ui/core/CircularProgress';
import SaveIcon from '@material-ui/icons/Save';
import Button from '@material-ui/core/Button';
import Typography from '@material-ui/core/Typography';
import Table from '@material-ui/core/Table';
import TableBody from '@material-ui/core/TableBody';
import TableCell from '@material-ui/core/TableCell';
import TableHead from '@material-ui/core/TableHead';
import TableRow from '@material-ui/core/TableRow';
import TableContainer from '@material-ui/core/TableContainer';
import Paper from '@material-ui/core/Paper';
import Chip from '@material-ui/core/Chip';
import FormControl from '@material-ui/core/FormControl';
import InputLabel from '@material-ui/core/InputLabel';
import Fuse from 'fuse.js';

const attributeStyles = makeStyles((theme) => ({
  interactionFields: {
    paddingTop: '30px'
  },
  tableContainer: {
    paddingTop: '30px'
  },
  tableLabel:{
    paddingTop: '18px',
    paddingLeft: '18px',
    width: 'max-content',
    userSelect: 'none'
  },
  table: {
    boxSizing: 'border-box',
    height: 'auto',
    overflowY: 'auto'
  },
  grayCell: {
    color: 'gray'
  }
}));

const fuseOptions = {
	isCaseSensitive: false,
	includeScore: true,
	shouldSort: true,
	includeMatches: false, // means char locations of match
	findAllMatches: true,
	minMatchCharLength: 2,
	keys: [
		'val'
	]
};

export type fuseSearchSets = {[OptionSet:string]: {val: string}[]};

type valueMatches = {[value: string]: {[OptionSet: string]: string[]}};

function SelectionTable({
  valueMatches,
  optionSetNames,
  selectedMatches,
  onSelect,
}: {
  valueMatches: valueMatches;
  optionSetNames: string[];
  selectedMatches: string[];
  onSelect: (m: string) => void;
}) {

  const classes = attributeStyles();

  const userValues = Object.keys(valueMatches);

  return <TableContainer component={Paper} className={classes.table}>
    <Table stickyHeader size='small'>
      <TableHead>
        <TableRow>
          <TableCell
            className={classes.grayCell}
            align='left'
            title='Query Value'
            key='header-query'
          >
            Query
          </TableCell>
          {optionSetNames.map((optionSetName) => <TableCell
              className={classes.grayCell}
              align='left'
              title={optionSetName}
              key={'header-optionSet-'+optionSetName}
            >
              {optionSetName}
            </TableCell>)
          }
        </TableRow>
      </TableHead>
      <TableBody>
        {userValues.length > 0 && (
          userValues.map((query, ind) => {
            return <TableRow key={'matches-row-'+ind}>
              <TableCell
                className={classes.grayCell}
                align='left'
                title={query}
                key={'matches-cell-query'}
              >
                {query}
              </TableCell>
              {optionSetNames.map((setName) => <TableCell
                  align='left'
                  title={query+'-'+setName+'-matches'}
                  key={'matches-cell-'+setName}
                >
                  {valueMatches[query][setName].map((match) => <Chip
                    key={query+'-'+setName+'-match-'+match}
                    size="small"
                    label={match}
                    color={selectedMatches.includes(match) ? 'primary' : 'default'}
                    onClick={() => onSelect(match)}
                  />)}
                </TableCell>
              )}
            </TableRow>;
          })
        )}
      </TableBody>
    </Table>
  </TableContainer>;
};

export default function SelectionsFromPasteModal({
  options,
  onSave,
  onClose,
  genomicData = true
}: {
  options: fuseSearchSets;
  onSave: (matches: string[]) => void;
  onClose: () => void;
  genomicData?: boolean;
}) {

  const classes = attributeStyles();

  // const classes = useStyles();
  const [userText, setUserText] = useState('');
  const [userValues, setUserValues] = useState([] as string[]);
  const [valueMatches, setValueMatches] = useState({} as valueMatches);
  const [selectedMatches, setSelectedMatches] = useState([] as string[]);
  const [busy, setBusy] = useState(false);

  const valuesNotMatched = useMemo(()=>{
    if (userValues.length==0) return [];
    return userValues.filter((val) => {
      const nearests = ([] as string[]).concat(...Object.values((valueMatches as valueMatches)[val]));
      return !nearests.some(nearest => selectedMatches.includes(nearest));
    });
  }, [userValues, valueMatches, selectedMatches]);

  // onClick for matches
  function addOrRemove(match: string) {
    if (selectedMatches.includes(match)) {
      setSelectedMatches(selectedMatches.filter(v => v!=match));
    } else {
      setSelectedMatches(selectedMatches.concat([match]));
    };
  };

  // onClick for Search button
  function parseAndMatchText() {
    setBusy(true);

    // Parse Input Text to discrete search elements
    const parsedValues = userText.split(/[\s\n,;]+/).filter((val) => !['',' '].includes(val));
    if (parsedValues.length<1) {
      setBusy(false);
      return;
    };
    setUserValues(parsedValues);

    // Perform fuzzy search per each options set
    const setMatches: {[OptionSet: string]: ({
      userValue: string,
      useMatch: boolean,
      bestMatch: string,
      topMatches: string[]
    } | {
      userValue: string,
      useMatch: false,
      bestMatch: null,
      topMatches: string[]
    })[]} = {};
    for (const set of Object.keys(options)) {
      const fuse = new Fuse(options[set], fuseOptions);

      setMatches[set] = parsedValues.map((val: string) => {
        let topMatches = fuse.search(val).slice(0,20);
        // Case: No matches picked at all
        if (topMatches.length<1) {
          return {
            userValue: val,
            useMatch: false,
            bestMatch: null,
            topMatches: [] as string[]
          };
        };
        // Case: No near-exact match.
        // For genomics data, pull thought-to-be desired top hits forward
        // Targets:
        //   - because of CD#s: letters, not numbers, following the query
        //   - because of how scData tools make markers unique: '.#' following the query.
        if (genomicData && topMatches[0].score as number >= 0.001) {
          const likely: number[] = [];
          const lesslikely: number[] = [];
          const valNoSymbol: string = val.replace(/\W/g,'.'); // to protect against symbols in the user input
          for (let ind = 0; ind < topMatches.length; ind++) {
            if (topMatches[ind].item.val.search(RegExp(valNoSymbol+'[a-z]|'+valNoSymbol+'\\.\\d', 'i')) > -1) {
              likely.push(ind);
            } else {
              lesslikely.push(ind);
            };
          };
          topMatches = likely.concat(lesslikely).map(i => topMatches[i]);
        };

        // Return (up to) top 5 matches
        return {
          userValue: val,
          useMatch: topMatches[0].score as number < 0.001,
          bestMatch: topMatches[0].item.val as string,
          topMatches: topMatches.map((v: any) => v.item.val as string).slice(0,5) as string[]
        };
      });
    };

    // Reformat output per userValue
    let matches = [] as string[];
    let valMatches: valueMatches = {};
    for (let ind = 0; ind < parsedValues.length; ind++) {
      const parsedValue = parsedValues[ind];
      const newMatches: {[OptionSet: string]: string[]} = {};
      for (const set of Object.keys(options)) {
        newMatches[set] = setMatches[set][ind].topMatches;
        if (setMatches[set][ind].useMatch) {
          if (!matches.includes(setMatches[set][ind].bestMatch as string)) {
            matches.push(setMatches[set][ind].bestMatch as string);
          };
        };
      };
      valMatches[parsedValue] = newMatches;
    };

    setValueMatches(valMatches);
    setSelectedMatches(matches);
    setBusy(false);
  }

  let howItWorksText: string= 'Input text is parsed into query values by splitting at spaces, newlines, commas, and semi-colons. Then each query value is compared to options of each Option Set of the data.';
  if (genomicData) {
    howItWorksText+=' Matches followed by a letter or ".#" are then prioritized over other near-matches to optimize the system for CD#s and how common genomic data tools deal with marker name overlaps.';
  };
  howItWorksText+=' Up to 5 top-ranked matches are then shown, per query, per Option Set.';

  return (
    <>
      <DialogTitle>
        Multi-Marker Add
      </DialogTitle>
      <DialogContent>
        <Typography key={'text1'}>
          <strong>{'Instructions: '}</strong>
          Paste or type into the Input box, then click the Search button. Results will then be shown in a table format where (near) perfect matches will be auto-selected.  Click any option to toggle between selected (orange) or deselected (gray).
        </Typography>
        <Typography key={'text2'}>
          <strong>{'How it works: '}</strong>
          {howItWorksText}
        </Typography>
        <Grid container spacing={1} className={classes.interactionFields}>
          <Grid item xs={3}>
            <TextField
              key={'TextField-input'}
              label={'Input'}
              InputLabelProps={{ shrink: true }}
              inputProps={{ 'data-testid': 'bulk-add-user-text-input' }}
              value={userText}
              onChange={(event: any) => {setUserText(event.target.value);}}
              autoFocus
              multiline
              size="small"
              variant='outlined'
              fullWidth
              placeholder={'Paste or type here'}
            />
          </Grid>
          <Grid item>
            <Button
              onClick={() => {parseAndMatchText();}}
              startIcon={busy ? <CircularProgress/> : <SearchIcon/>}
              disabled={busy}
              color='secondary'
              variant='contained'
            >
              Search
            </Button>
          </Grid>
          {Object.keys(valueMatches).length < 1 ? 
            <Grid item xs={7} container direction='column'>
              <TextField
                label={'Output Area'}
                InputLabelProps={{ shrink: true }}
                multiline
                value={''}
                size="small"
                variant='outlined'
                disabled
                fullWidth
                placeholder={'Awaiting Search...'}
              />
            </Grid> :
            <Grid item xs={7} container direction='column'>
              <Grid item>
                <TextField
                  key={'TextField-output'}
                  label={'Values not Matched'}
                  InputLabelProps={{ shrink: true }}
                  multiline
                  value={valuesNotMatched.join(', ')}
                  size="small"
                  disabled
                  error={valuesNotMatched.length!=0}
                  variant='outlined'
                  fullWidth
                  placeholder={'None!'}
                />
              </Grid>
              <Grid item>
                <FormControl className={classes.tableContainer}>
                  <InputLabel shrink
                    className={classes.tableLabel}
                  >
                    Value Match Options
                  </InputLabel>
                  <SelectionTable
                    valueMatches={valueMatches}
                    optionSetNames={Object.keys(options)}
                    selectedMatches={selectedMatches}
                    onSelect={(match: string)=>{addOrRemove(match);}}
                  />
                </FormControl>
              </Grid>
            </Grid>
          }
        </Grid>
      </DialogContent>
      <DialogActions>
        <Button
          onClick={() => {onSave(selectedMatches); onClose();}}
          startIcon={<SaveIcon />}
          disabled={selectedMatches==null}
          color='primary'
          variant='contained'
          aria-label='Add Selected Matches'
        >
          Add Selected Matches
        </Button>
        <Button
          onClick={() => {onClose();}}
          startIcon={<CancelIcon />}
          color='secondary'
          variant='contained'
          aria-label='Cancel'
        >
          Cancel
        </Button>
      </DialogActions>
    </>
  );
};
