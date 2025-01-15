import React, {
  useState,
  useCallback,
  useContext,
  useEffect,
  useMemo
} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';

import {
  inputGroupName,
  missingConfigInputDefaults
} from '../../../selectors/workflow_selectors';
import PrimaryInputGroup from './primary_input_group';
import {
  bindInputSpecification,
  BoundInputSpecification,
  getParamUISpecifications,
  InputSpecification
} from '../user_interactions/inputs/input_types';
import {useWorkspace} from '../../../contexts/workspace_context';
import {Maybe, maybeOfNullable} from '../../../selectors/maybe';
import {
  BufferedInputsContext,
  WithBufferedInputs
} from '../../../contexts/input_state_management';
import Collapse from '@material-ui/core/Collapse';
import Typography from '@material-ui/core/Typography';
import Grid from '@material-ui/core/Grid';
import Card from '@material-ui/core/Card';
import CardContent from '@material-ui/core/CardContent';
import ExpandMoreIcon from '@material-ui/icons/ExpandMore';
import ExpandLessIcon from '@material-ui/icons/ExpandLess';
import IconButton from '@material-ui/core/IconButton';
import {makeStyles} from '@material-ui/core/styles';

const useStyles = makeStyles((theme) => ({
  card: {
    borderRadius: 0,
    border: '1px solid #eee'
  },
  header: {
    padding: '5px 10px',
    cursor: 'pointer'
  }
}));

export default function PrimaryInputs() {
  const {state} = useContext(VulcanContext);
  const {workspace} = state;
  if (!workspace) return null;

  const inputSpecifications = useMemo(
    () => getParamUISpecifications(workspace),
    [workspace]
  );

  let groupedInputs = useMemo(() => {
    return inputSpecifications.reduce((result, spec) => {
      let groupName = inputGroupName(spec.name) || 'Config Inputs';
      result[groupName] = result[groupName] || [];
      result[groupName].push(spec);
      return result;
    }, {} as {[k: string]: InputSpecification[]});
  }, [
    inputSpecifications
  ]);

  const [expanded, setExpanded] = useState([] as string[]);
  function toggleExpanded(groupName: string) {
    if (expanded.includes(groupName)) {
      setExpanded(expanded.filter(val => val != groupName));
    } else {
      setExpanded(expanded.concat(groupName))
    }
  }

  const [open, setOpen] = useState(false);

  const classes = useStyles();

  return (
    <Card className={classes.card}>
      <Grid
        className={classes.header}
        container
        alignItems='center'
        justifyContent='space-between'
        onClick={() => setOpen(!open)}
      >
        <Grid item>
          <Typography variant='h6'>Config Params</Typography>
        </Grid>
        <Grid item>
          <IconButton size='small'>
            {open ? (
              <ExpandLessIcon fontSize='small' />
            ) : (
              <ExpandMoreIcon fontSize='small' />
            )}
          </IconButton>
        </Grid>
      </Grid>
      <Collapse in={open}>
        {Object.keys(groupedInputs)
          .sort()
          .map((groupName, index) => {
            return (
              <PrimaryInputGroup
                expanded={expanded.includes(groupName)}
                select={() => toggleExpanded(groupName)}
                groupName={groupName.split('_').slice(1).join(' ')}
                key={index}
                inputs={groupedInputs[groupName]}
              />
            );
          })}
      </Collapse>
    </Card>
  );
}
