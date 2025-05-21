import React, {useMemo} from 'react';
import {DataEnvelope, WithInputParams} from '../../input_types';
import {maybeOfNullable} from '../../../../../selectors/maybe';
import {useSetsDefault} from '../../useSetsDefault';
import { SingleMagmaRecordSelectionPieceRct } from '../pieces/magma_pickers';

declare const CONFIG: {[key: string]: any};

function getFromData(key: string, data: DataEnvelope<any>, showError: (e:string) => void) {
  return useMemo(()=>{
    if (key in data) {
      return data[key]
    } else {
      showError(`Target ${key} not provided`)
      return undefined
    }
  }, [data])
}

export function MagmaRecordInput({onChange, label, defaultValue, showError, data, ...props}: WithInputParams<{}, string, string[] >
) {
  const value = useSetsDefault(defaultValue || null, props.value, onChange, 'record');
  const modelName = getFromData('modelName', data, (e) => showError(e))
  const targetAttribute = getFromData('targetAtt', data, (e) => showError(e))
  const showAttributeNames = getFromData('otherAttsShow', data, (e) => showError(e))
  const hasAttributeNames = getFromData('hasAtts', data, (e)=>{}) // Optional

  return <SingleMagmaRecordSelectionPieceRct
      name={`magma-record-selection-${label}`}
      changeFxn={(v, k) => {onChange({record: maybeOfNullable(v)})}}
      value={value}
      label={label || `${modelName} Record`}
      modelName={modelName}
      targetAttribute={targetAttribute}
      hasAttributes={hasAttributeNames}
      attributesDisplay={showAttributeNames}
      showError={showError}
      projectName={CONFIG.project_name}
    />
};