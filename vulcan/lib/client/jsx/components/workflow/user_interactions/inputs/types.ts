import React from 'react';

export type InputType = string;

export interface InputSpecification {
    type: InputType,
    label: string | undefined,
    name: string,
    default: any | null,
    data?: {[k: string]: any} | undefined,
    doc?: string | undefined,
}

export type InputOnChange = (inputName: string, val: any) => void

export type InputBackendComponent = (p: { input: InputSpecification, onChange: InputOnChange}) => React.ReactElement | null;
