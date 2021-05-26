import React, {createContext} from 'react';
import {Workflow} from "../api_types";

const defaultWorkflowContext = {} as Workflow;
export const WorkflowContext = createContext(defaultWorkflowContext);
