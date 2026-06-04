import React, {
  useCallback,
  useContext,
  useEffect,
  useMemo,
  useState
} from 'react';
import 'regenerator-runtime/runtime';

import Grid from '@material-ui/core/Grid';

import {makeStyles} from '@material-ui/core/styles';

import {VulcanContext} from '../../contexts/vulcan_context';
import {Workflow} from '../../api_types';
import WorkflowCard from './workflow_card';
const useStyles = makeStyles((theme) => ({
  workflows: {
    width: '100%',
    padding: '5px',
    boxShadow: '0 0 15px 0px #f5f5f5 inset'
  }
}));

export default function WorkflowsCarousel({
  project_name,
  onSelectWorkflow
}: {
  project_name: string;
  onSelectWorkflow: (workflow: Workflow | null) => void;
}) {
  const [selectedWorkflow, setSelectedWorkflow] = useState<Workflow | null>(
    null
  );

  const classes = useStyles();
  let {state} = useContext(VulcanContext);
  const {workflows} = state;

  // v2: workflows are already retrieved in a project specific context
  const projectWorkflows = !!workflows
    ? workflows
    : [];

  const collator = new Intl.Collator(undefined, {
    numeric: true,
    sensitivity: 'base'
  });
  const sortedWorkflows = useMemo(() => {
    return projectWorkflows.sort((a, b) =>
      collator.compare(a.name, b.name)
    );
  }, [projectWorkflows]);

  useEffect(() => {
    onSelectWorkflow(selectedWorkflow);
  }, [selectedWorkflow]);

  const handleOnSelectWorkflow = useCallback(
    (workflow: Workflow) => {
      if (workflow.id === selectedWorkflow?.id) {
        setSelectedWorkflow(null);
      } else {
        setSelectedWorkflow(workflow);
      }
    },
    [selectedWorkflow]
  );

  return (
    <Grid
      className={classes.workflows}
      container
    >
      {sortedWorkflows.map((workflow: Workflow, index: number) => {
        return (
          <Grid item key={index}>
            <WorkflowCard
              selected={selectedWorkflow?.id === workflow.id}
              workflow={workflow}
              onClick={() => handleOnSelectWorkflow(workflow)}
            />
          </Grid>
        );
      })}
    </Grid>
  );
}
