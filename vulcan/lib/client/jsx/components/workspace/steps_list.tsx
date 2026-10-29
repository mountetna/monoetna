import React, {useContext, useState, useEffect, useMemo} from 'react';

import {makeStyles} from '@material-ui/core/styles';
import {VulcanContext} from '../../contexts/vulcan_context';
import {outputUINames, outputUIsWithInputsReady, stepOfName} from '../../selectors/workflow_selectors';
import {useWorkspace} from '../../contexts/workspace_context';
import { VulcanConfig } from '../../api_types';
import StepIconName from './drawers/step_elements/step_icon_name';
import FlatButton from 'etna-js/components/flat-button'
import Grid from '@material-ui/core/Grid';
import Drawer from '@material-ui/core/Drawer';
import IconButton from '@material-ui/core/IconButton';
import List from '@material-ui/core/List';
import ListItem from '@material-ui/core/ListItem';
import Typography from '@material-ui/core/Typography';
import ChevronRightIcon from '@material-ui/icons/ChevronRight';
import ChevronLeftIcon from '@material-ui/icons/ChevronLeft';

const useStyles = makeStyles((theme) => ({
  steps: {
    borderLeft: '1px solid #ccc'
  },
  stepsOpen: {
    width: '800px',
    overflow: 'auto'
  },
  stepsClosed: {
    width: '60px'
  },
  stepRow: {
    display: 'flex',
    gap: '5px'
  },
  step: {
    border: '1px solid red',
    height: '60px',
    display: 'flex',
    padding: '5px 10px',
    background: '#f2eeee',
    alignItems: 'center',
    borderRadius: '2px',
    position: 'relative',
    cursor: 'pointer',
    boxShadow: '0 0 2px #ccc',
    '&:hover': {
      background: '#fee',
    }
  },
  stepLabel: {
    width: '200px',
    overflow: 'hidden',
    display: 'block',
    textOverflow: 'ellipsis',
    textWrap: 'nowrap'
  },
  connections: {
    position: 'absolute',
    display: 'flex',
    gap: 5,
    right: 0,
    width: 'calc(100% - 10px)',
    justifyContent: 'flex-end',
    padding: '5px'
  },
  input: {
    width: '10px',
    height: '10px',
    position: 'relative',
    boxShadow: '0 0 2px #444',
    '&:hover': {
      '& $wire path': {
        strokeWidth: 5,
        strokeOpacity: 1
      }
    }
  },
  resource: {
    background: '#88f',
    borderRadius: '50%',
  },
  pending: {
    background: 'red',
  },
  missing: {
    width: '8px',
    height: '8px',
    border: '1px solid red',
    background: 'none',
  },
  ok: {
    background: '#1ce31c',
    borderRadius: '50%'
  },
  output: {
    width: '10px',
    height: '10px',
    boxShadow: '0 0 2px #444',
  },
  wire: {
    position: 'absolute',
    bottom: 0,
    pointerEvents: 'none'
  }
}));

const snakemakeSummaryStyle = (summary, classes) => {
  if (!summary) return classes.resource;
  if (summary.status == "ok") return classes.ok;
  if (summary.status == "missing") {
    if (summary.plan == "update pending") return classes.pending;
    return classes.missing;
  }
  return classes.pending;
}

const FileOutput = ({output, summary}) => {
  const classes = useStyles();
  let { status, plan } = summary || { status: 'resource' }
  let outputStyle = snakemakeSummaryStyle(summary, classes);
  return <div className={ `${outputStyle} ${classes.output}` } title={ `${output} - ${status}${ plan ? ` - ${plan}` : '' }` }/>
}

const FileInput = ({input, pos, stepPos, stepRow, wireTo, summary}) => {
  const classes = useStyles();
  const PAD = 5;

  let inputStyle = snakemakeSummaryStyle(summary, classes);
  let { status, plan } = summary || { status: 'resource' }

  const [ path, bounds ] = wireTo( pos, stepPos, stepRow, input );
  return <>
    <div className={ `${inputStyle} ${classes.input}` } title={ `${input} - ${status}${ plan ? ` - ${plan}` : '' }` }>
    {
      bounds && <div className={ classes.wire } style={{ right: bounds.flip ? undefined : 0, left: bounds.flip ? 0 : undefined }}>
      <svg width={ bounds.width + 2*PAD } height={ bounds.height + 2*PAD } viewBox={ `${bounds.x-PAD} ${bounds.y-PAD} ${bounds.width+2*PAD} ${bounds.height+2*PAD}` }>
        <path fill="none" stroke="goldenrod" strokeOpacity="0.1" strokeWidth="2" d={ path }/>
      </svg>
    </div>
    }
    </div>
  </>
}



const Step = ({name, files, summary, config, wireTo, stepPos, row}) => {
  const classes = useStyles();

  const label = config ? config.display : name.replace(/_/g,' ');

  return <div className={ classes.step }>
    <div className={ classes.stepLabel } title={label}>{label}</div>
    {
      files.outputs.length > 0 && <div className={ classes.connections} style={{ bottom: 0 }}>
        {
          files.outputs.map(
            output => <FileOutput key={output} summary={summary[output]} output={ output }/>
          )
        }
      </div>
    }
    {
      files.inputs.length > 0 && <div className={ classes.connections} style={{ top: 0, background: "#fda" }}>
        {
          files.inputs.map(
            (input,pos) => <FileInput key={input} 
              input={ input }
              pos={files.inputs.length - pos - 1}
              summary={summary[input]}
              stepPos={stepPos}
              stepRow={row}
              wireTo={wireTo}/>
          )
        }
      </div>
    }
  </div>
}

export default function StepsList() {
  const [open, setOpen] = useState(true);
  const {state} = useContext(VulcanContext);
  const {status, workspace, refreshingProgress} = state;

  const classes = useStyles();

  function handleToggle() {
    setOpen(!open);
  }

  if (!workspace || Object.keys(workspace).length == 0) return null;

  let rows = []
  let rowPos = {}
  const visited = new Set([]);
  let currentRow = 0;

  const { dag } = workspace;

  while ( (visited.size < Object.keys(dag).length)) {
    rows[currentRow] = []

    const justVisited = [];
    Object.keys(dag).forEach(
      step => {
        if (!visited.has(step) && dag[step].every( s => visited.has(s)))
          {
            rows[currentRow].push(step);
            workspace.filegraph[step].outputs.forEach(
              (output,i) => {
                rowPos[output] = [
                  currentRow,
                  rows[currentRow].length - 1,
                  workspace.filegraph[step].outputs.length - i - 1
                ];
              }
            );
            justVisited.push(step);
          }
      }
    )
    for (let step of justVisited) visited.add(step)
    currentRow += 1;
  }

  console.log({rowPos});

  const wireTo = ( fromPos, fromStepPos, fromRow, input ) => {
    if (!rowPos[input]) return [ '', null ];
    const [ toRow, toStepPos, toPos ] = rowPos[input];

    const NODE_GAP = 15;
    const yDiff = (toRow - fromRow) * 86 + 57;
    const xDiff = (toStepPos - fromStepPos) * 227 - (toPos - fromPos) * NODE_GAP;
    const startX = (xDiff > 0) ? -6 : 0;
    const startY = 0;

    return [
      `M ${startX},${startY} c ${0},${yDiff}, ${xDiff},${0} ${xDiff},${yDiff}`,
      { width: Math.abs(xDiff) - startX, height: Math.abs(yDiff) - startY, x: (xDiff > 0) ? startX : startX + xDiff, y: startY + yDiff, flip: xDiff > 0 }
    ]
  }

  return (
    <Grid className={ `${classes.steps} ${open ? classes.stepsOpen : classes.stepsClosed}` } >
      <Grid alignItems='center' container>
        <IconButton onClick={ () => setOpen(!open) }>
          { open ? <ChevronRightIcon/> : <ChevronLeftIcon/> }
        </IconButton>
      </Grid>
      <List>
        {rows.map(
          (row,i) => <ListItem key={i}>
            <div className={ classes.stepRow }>
            {
              row.map( (stepName,pos) => <Step key={stepName} row={i} stepPos={pos} name={stepName} files={workspace.filegraph[stepName]} summary={workspace.summary} config={workspace.vulcan_config[stepName]} wireTo={wireTo}/>
              )
            }
            </div>
          </ListItem>
        )}
      </List>
    </Grid>
  );
};
