import React, {useContext, useState, useEffect, useMemo} from 'react';

import {makeStyles} from '@material-ui/core/styles';
import {VulcanContext} from '../../contexts/vulcan_context';
import {outputUINames, outputUIsWithInputsReady, stepOfName} from '../../selectors/workflow_selectors';
import {useWorkspace} from '../../contexts/workspace_context';
import { VulcanConfig } from '../../api_types';
import StepIcon from './drawers/step_elements/step_icon';
import FlatButton from 'etna-js/components/flat-button'
import Grid from '@material-ui/core/Grid';
import Drawer from '@material-ui/core/Drawer';
import IconButton from '@material-ui/core/IconButton';
import Tooltip from '@material-ui/core/Tooltip';
import List from '@material-ui/core/List';
import ListItem from '@material-ui/core/ListItem';
import Typography from '@material-ui/core/Typography';
import ChevronRightIcon from '@material-ui/icons/ChevronRight';
import ChevronLeftIcon from '@material-ui/icons/ChevronLeft';
import LiveHelpIcon from '@material-ui/icons/LiveHelp';
import ExploreIcon from '@material-ui/icons/Explore';
import CloseIcon from '@material-ui/icons/Close';
import SettingsApplicationsIcon from '@material-ui/icons/SettingsApplications';

const useStyles = makeStyles((theme) => ({
  steps: {
    borderLeft: '1px solid #ccc',
    height: '100%',
    overflow: 'hidden'
  },
  stepsHeader: {
    borderBottom: '1px solid #eee'
  },
  stepsOpen: {
    width: '800px',
  },
  stepsClosed: {
    width: '48px'
  },
  stepsLegendIcon: {
    color: '#888', position: 'absolute', top: '95px', right: '10px',
    cursor: 'pointer',
    '&:hover': {
      color: '#ccc'
    }
  },
  stepsLegendDismiss: {
    color: '#888', position: 'absolute', top: '5px', right: '5px',
    cursor: 'pointer',
    '&:hover': {
      color: '#ccc'
    }
  },
  stepsLegend: {
    position: 'absolute', top: '100px', right: '10px',
    border: '1px solid #888',
    boxShadow: '0 0 2px #ccc',
    display: 'flex',
    flexDirection: 'column',
    zIndex: 2,
    background: 'rgb(from white r g b / 85%)',
    width: '110px',
    padding: '10px'
  },
  stepGraph: {
    overflow: 'auto',
    padding: '15px',
    display: 'flex',
    flexDirection: 'column',
    height: 'calc( 100% - 49px - 30px )',
    gap: '25px'
  },
  stepList: {
    overflow: 'auto',
    height: '100%'
  },
  stepRow: {
    display: 'flex',
    gap: '25px'
  },
  stepSummary: {
    width: '48px',
    height: '47px',
    borderBottom: '1px solid #eee',
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'center'
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
    }
  },
  stepLabel: {
    width: '200px',
    overflow: 'hidden',
    display: 'flex',
    textWrap: 'nowrap',
    position: 'relative',
    gap: '5px',
    '&:before': {
      content:"''",
      width: '100%',
      height: '100%',
      position: 'absolute',
      pointerEvents: 'none',
      left:0,
      top:0,
      background: 'linear-gradient(to right, transparent 180px, #f2eeee)'
    }
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
  infobox: {
    position: 'absolute',
    display: 'flex',
    gap: 5,
    left: 0,
    top: 0,
    color: '#888',
    justifyContent: 'flex-start'
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
  legendItem: {
    display: 'flex',
    alignItems: 'center',
    gap: '5px',
    '& div': {
      width: '10px',
      height: '10px'
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
  feed: {
    borderRadius: '2px',
    border: '1px solid purple',
    transform: 'rotate(45deg)'
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

const FileOutput = ({output, summary, highlight, setHighlight, feedOutput}) => {
  const classes = useStyles();
  let { status, plan } = summary || { status: 'resource' }
  let outputStyle = snakemakeSummaryStyle(summary, classes);
  return <Tooltip title={ `${output} - ${status}${ plan ? ` - ${plan}` : '' }` }>
    <div className={ `${outputStyle} ${classes.output} ${feedOutput ? classes.feed : ''}` } onMouseOver={ () => setHighlight(new Set([output]))} onMouseOut={ () => setHighlight(new Set([])) }/>
  </Tooltip>
}

const FileInput = ({input, pos, stepPos, stepRow, wireTo, summary, highlight}) => {
  const classes = useStyles();
  const PAD = 5;

  let inputStyle = snakemakeSummaryStyle(summary, classes);
  let { status, plan } = summary || { status: 'resource' }

  const [ path, bounds ] = wireTo( pos, stepPos, stepRow, input );
  return <Tooltip title={ `${input} - ${status}${ plan ? ` - ${plan}` : '' }` }>
    <div className={ `${inputStyle} ${classes.input}` }>
    {
      bounds && <div className={ classes.wire } style={{ right: bounds.flip ? undefined : 0, left: bounds.flip ? 0 : undefined }}>
        <svg width={ bounds.width + 2*PAD } height={ bounds.height + 2*PAD } viewBox={ `${bounds.x-PAD} ${bounds.y-PAD} ${bounds.width+2*PAD} ${bounds.height+2*PAD}` }>
          <path fill="none" stroke="goldenrod" strokeOpacity={ highlight.has(input) ? 1 : 0.1 } strokeWidth={ highlight.has(input) ? 5 : 2 } d={ path }/>
        </svg>
      </div>
    }
    </div>
  </Tooltip>
}

const UIStepIndicator = () => <Tooltip title="This step requires user input"><SettingsApplicationsIcon fontSize="small"/></Tooltip>

const StepsLegend = () => {
  const [ shown, setShown ] = React.useState(false);
  const classes = useStyles();
  if (!shown) return <div className={ classes.stepsLegendIcon } onClick={ () => setShown(true) }>
    <ExploreIcon/>
  </div>

  return <div className={ `${classes.stepsLegend}` }>
    <IconButton size="small" onClick={ () => setShown(false) } className={ classes.stepsLegendDismiss }>
      <CloseIcon fontSize="small"/>
    </IconButton>
    <div className={ classes.legendItem }><div className={ classes.missing}/> missing</div>
    <div className={ classes.legendItem }><div className={ classes.pending}/> pending</div>
    <div className={ classes.legendItem }><div className={ classes.ok}/> present</div>
    <div className={ classes.legendItem }><div className={ classes.resource}/> resource</div>
    <div className={ classes.legendItem }><div className={ `${classes.feed} ${classes.ok}` }/> output</div>
  </div>
}

const StepSummary = ({name, files, summary, config}) => {
  const classes = useStyles();

  return <div className={ classes.stepSummary }>
    <StepIcon step={
      stepOfName(name, config ? { [name]: config } : {} )
    } showName={true}/>
  </div>
}

const Step = ({name, files, summary, config, wireTo, stepPos, row, highlight, setHighlight,outputFeedFiles}) => {
  const classes = useStyles();

  const label = config ? config.display : name.replace(/_/g,' ');

  return <div className={ classes.step }>
    <div className={ classes.stepLabel }>
      <StepIcon step={
        stepOfName(name, config ? { [name]: config } : {} )
      }/>
      <Tooltip title={label}><span>{label}</span></Tooltip>
    </div>
    {
      files.outputs.length > 0 && <div className={ classes.connections} style={{ bottom: 0 }}>
        {
          files.outputs.map(
            output => <FileOutput key={output} summary={summary[output]} output={ output } highlight={highlight} feedOutput={ outputFeedFiles.has(output) } setHighlight={setHighlight}/>
          )
        }
      </div>
    }
    {
      files.inputs.length > 0 && <div className={ classes.connections} style={{ top: 0, background: "#fda" }}>
        {
          files.inputs.map(
            (input,pos) => <FileInput key={input} highlight={highlight} setHighlight={setHighlight} 
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
    <div className={ classes.infobox }>
      { config?.name ? <UIStepIndicator/> : null }
    </div>
  </div>
}

export default function StepsList() {
  const [open, setOpen] = useState(false);
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

  const outputFeedFiles = new Set([]);

  Object.keys(workspace.vulcan_config).forEach(
    outputName => {
      const config = workspace.vulcan_config[outputName];
      if (config.output) return null;
      if (!config?.input?.files) return null;
      config.input.files.forEach(
        file => outputFeedFiles.add(`output/${file}`)
      )
    }
  );

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
    const yDiff = (toRow - fromRow) * (72 + 25) + 59;
    const xDiff = (toStepPos - fromStepPos) * (222 + 25) - (toPos - fromPos) * NODE_GAP;
    const startX = (xDiff > 0) ? -6 : 0;
    const startY = 0;

    return [
      `M ${startX},${startY} c ${0},${yDiff}, ${xDiff},${0} ${xDiff},${yDiff}`,
      { width: Math.abs(xDiff) - startX, height: Math.abs(yDiff) - startY, x: (xDiff > 0) ? startX : startX + xDiff, y: startY + yDiff, flip: xDiff > 0 }
    ]
  }

  const [ highlight, setHighlight ] = useState(new Set([]));

  return (
    <Grid className={ `${classes.steps} ${open ? classes.stepsOpen : classes.stepsClosed}` } >
      { open && <StepsLegend/> }
      <Grid alignItems='center' container className={ classes.stepsHeader }>
        <IconButton onClick={ () => setOpen(!open) }>
          { open ? <ChevronRightIcon/> : <ChevronLeftIcon/> }
        </IconButton>
      </Grid>
      {
        open ?
        <div className={ classes.stepGraph }>
          {
            rows.map(
              (row,i) => <div key={i} className={ classes.stepRow }>
                {
                  row.map(
                    (stepName,pos) => <Step
                      key={stepName}
                      highlight={highlight}
                      setHighlight={setHighlight}
                      row={i}
                      stepPos={pos}
                      name={stepName}
                      files={workspace.filegraph[stepName]}
                      summary={workspace.summary}
                      config={workspace.vulcan_config[stepName]}
                      outputFeedFiles={outputFeedFiles}
                      wireTo={wireTo}/>
                  )
                }
                </div>
            )
          }
        </div> : <div className={ classes.stepList }>
          {
            workspace.dag_flattened.map(
            (stepName,pos) => <StepSummary key={ stepName } name={stepName} files={workspace.filegraph[stepName]} summary={workspace.summary} config={workspace.vulcan_config[stepName]}/>
            )
          }
        </div>
      }
    </Grid>
  );
};
