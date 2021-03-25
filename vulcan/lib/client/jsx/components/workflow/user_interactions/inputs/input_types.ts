import * as React from 'react';
import {WorkflowStep} from "../../../../api_types";

export type InputType = string;

export interface InputSpecification {
    type: InputType,
    label: string | undefined,
    name: string, // output source
    default: any | null,
    data?: {[k: string]: any} | undefined,
    doc?: string | undefined,
}

export interface UIStep {
    step: WorkflowStep | GroupedInputStep,
    // Index in the original workflow
    index: number,
}

export type GroupedInputStep = WorkflowStep & {
    isGroup: true
}

export type InputOnChange = (inputName: string, val: any) => void

export type InputBackendComponent = (p: { input: InputSpecification, onChange: InputOnChange}) => React.ReactElement | null;
