import React from 'react';

export type InputType = 'int' | 'float' | 'boolean' | 'string' | 'array' | 'File' | 'multiselect-string' | 'select-autocomplete' | 'checkboxes';

export interface InputSpecification {
    type: InputType,
    label: string | undefined,
    name: string,
    default: any | null,
    data: {[k: string]: any},
    doc: string | undefined,
}

export type InputOnChange = (inputName: string, val: any) => void

export type InputBackendComponent = (p: { input: InputSpecification, onChange: InputOnChange}) => React.ReactElement | null;
