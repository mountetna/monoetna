import React, {useCallback, useContext, useState} from 'react';
import 'regenerator-runtime/runtime';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {pushLocation} from 'etna-js/actions/location_actions';

import {VulcanContext} from '../contexts/vulcan_context';
import Card from '../components/dashboard/card';
import {workflowName} from "../selectors/workflow_selectors";
import ReactModal from "react-modal";
import SelectInput from 'etna-js/components/inputs/select_input'

const modalStyles = {
  content: {
    top: '50%',
    left: '50%',
    right: 'auto',
    bottom: 'auto',
    marginRight: '-50%',
    transform: 'translate(-50%, -50%)',
    width: 400,
  }
};

export default function Dashboard() {
  const invoke = useActionInvoker();
  let {state} = useContext(VulcanContext);
  const {workflows} = state;
  const [selectedWorkflow, setSelectedWorkflow] = useState(null);
  const [selectedProject, setSelectedProject] = useState(null);
  const [projectModalOpen, setProjectModalOpen] = useState(false);

  const visitWorkflow = useCallback(() => {
    if (!selectedProject || !selectedWorkflow) return;
    invoke(pushLocation(`/workflow/${selectedProject}/${workflowName(selectedWorkflow)}`));
  }, [selectedWorkflow, selectedProject])

  const openProjectModal = useCallback((workflow) => {
    setSelectedWorkflow(workflow);
    setProjectModalOpen(true);
  }, []);

  const confirmProject = useCallback(() => {
    setProjectModalOpen(false);
    visitWorkflow();
  }, [visitWorkflow]);

  const cancelProjectModal = useCallback(() => {
    setSelectedWorkflow(null);
    setSelectedProject(null);
    setProjectModalOpen(false);
  }, []);

  if (!workflows || workflows.length === 0) return null;

  return (
    <main className='vulcan-dashboard'>
      {workflows.map((w, ind) => {
        return (
          <Card
            workflow={w}
            key={ind}
            onClick={() => {
              openProjectModal(w);
            }}
          />
        );
      })}
      { selectedWorkflow && <ReactModal
        isOpen={projectModalOpen}
        onRequestClose={cancelProjectModal}
        style={modalStyles}
        contentLabel='Select Your Project'
      >
        <div className="">
          To begin working on the {workflowName(selectedWorkflow)} workflow, you must select a project
          context to begin your work in.
        </div>
        <div className="select-project-modal">
          <SelectInput
            showNone
            defaultValue={null}
            values={selectedWorkflow.projects}
            onChange={setSelectedProject}
          />

          <button className='modal-button' onClick={confirmProject}>
            Confirm
          </button>
        </div>
      </ReactModal> }
    </main>
  );
}
