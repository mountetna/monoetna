import React, {useState, useCallback, useContext, useEffect, useMemo} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';

import {inputGroupName} from '../../../selectors/workflow_selectors';
import InputGroup from './input_group';
import {
  bindInputSpecification, BoundInputSpecification, DataEnvelope, getInputSpecifications
} from "../user_interactions/inputs/input_types";
import {useWorkflow} from "../../../contexts/workflow_context";
import {Maybe, maybeOfNullable} from "../../../selectors/maybe";
import {BufferedInputsContext, WithBufferedInputs} from "../../../contexts/input_state_management";
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
  const {commitSessionInputChanges, dispatch} = useContext(VulcanContext);

  return <WithBufferedInputs commitSessionInputChanges={commitSessionInputChanges} dispatch={dispatch} stepName={null}>
    <PrimaryInputsInner/>
  </WithBufferedInputs>
}


function PrimaryInputsInner() {
  const {state} = useContext(VulcanContext);
  const {inputs, setInputs} = useContext(BufferedInputsContext);
  const {session} = state;
  const {workflow} = useWorkflow();

  // Ensure defaults are set.
  useEffect(() => {
    let withDefaults: DataEnvelope<Maybe<any>> = {};
    Object.keys(workflow.inputs).forEach(inputName => {
      if (!(inputName in session.inputs) && !(inputName in inputs)) {
        withDefaults[inputName] = maybeOfNullable(workflow.inputs[inputName].default);
      }
    })

    if (Object.keys(withDefaults).length > 0) {
      console.log('setting buffered inputs', {...inputs, ...withDefaults})
      setInputs(inputs => ({...inputs, ...withDefaults}));
    }
  }, [inputs, session.inputs, setInputs, workflow.inputs])

  const inputSpecifications = useMemo(() =>
    getInputSpecifications(Object.entries(workflow.inputs), workflow), [workflow]);

  let groupedInputs = useMemo(() => {
    return inputSpecifications.reduce((result, spec) => {
      let groupName = inputGroupName(spec.name) || "Inputs";
      result[groupName] = result[groupName] || [];
      result[groupName].push(bindInputSpecification(spec,
        workflow,
        state.status,
        state.session,
        state.data,
        inputs,
        setInputs,
      ));
      return result;
    }, {} as {[k: string]: BoundInputSpecification[]});
  }, [inputSpecifications, inputs, setInputs, state.data, state.session, state.status, workflow]);

  console.log('in primary inputs', inputs, session, groupedInputs)
  const [ expanded, setExpanded ] = useState('');

  const [ open, setOpen ] = useState(false);

  const classes = useStyles();

  return (
    <Card className={classes.card}>
        <Grid className={classes.header} container alignItems='center' justify='space-between' onClick={ () => setOpen(!open) }>
          <Grid item>
            <Typography variant='h6'>Primary Inputs</Typography>
          </Grid>
          <Grid item>
            <IconButton size='small'>
              { open ? <ExpandLessIcon fontSize='small'/> : <ExpandMoreIcon fontSize='small'/>}
            </IconButton>
          </Grid>
        </Grid>
        <Collapse in={ open }>
          {Object.keys(groupedInputs)
            .sort()
            .map((groupName, index) => {
              return (
                <InputGroup
                  expanded={ expanded == groupName }
                  select={ () => setExpanded(expanded == groupName ? '' : groupName) }
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
