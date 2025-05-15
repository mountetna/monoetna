import React, { useState, useCallback, useEffect, useMemo } from 'react';
import {getModels, getDocuments} from 'etna-js/api/magma_api'
import Autocomplete from '@material-ui/lab/Autocomplete';
import { Grid, IconButton, TextField } from '@material-ui/core';
import { json_get } from 'etna-js/utils/fetch';
import {
  arrayLevels,
  DisabledTextbox,
  PieceBaseInputs
} from './user_input_pieces';
import { DropdownPieceRct } from './dropdown_piece';

declare const CONFIG: {[key: string]: any};

interface SingleMagmaRecordSelectionInputs extends PieceBaseInputs<string | null> {
  modelName: string;
  targetAttribute: string;
  hasAttributes?: string[];
  attributesDisplay?: string[];
  showError?: (e:string) => void;
  projectName: string;
};

export function SingleMagmaRecordSelectionPieceRct({
  name,
  changeFxn,
  value,
  label,
  modelName,
  targetAttribute,
  hasAttributes,
  attributesDisplay = [targetAttribute],
  showError = (e) => {},
  projectName
}: SingleMagmaRecordSelectionInputs): React.ReactElement | null {

  const [allRecords, setAllRecords] = useState({})
  const [filtRecords, setfiltRecords] = useState([] as string[])
  useEffect(() => {
    const targets = !hasAttributes ? 
      arrayLevels(attributesDisplay.concat(targetAttribute)) :
      arrayLevels(attributesDisplay.concat(hasAttributes).concat(targetAttribute));
    getDocuments({
      project_name: projectName,
      model_name: modelName,
      record_names: 'all',
      attribute_names: targets
    }).then((recs) => {
      const allRecs = recs.models[modelName].documents
      setAllRecords(recs.models[modelName].documents)
    }).catch((e) => showError(e));
  }, [modelName, attributesDisplay, hasAttributes])

  useEffect(() => {
    if (!!hasAttributes && Object.keys(allRecords).length>0) {
      const filt = Object.keys(allRecords).filter((v) => {
        return hasAttributes.every((has) => allRecords[v][has]!=null)
      })
      setfiltRecords(filt);
    } else {
      Object.keys(allRecords)
    }
  }, [allRecords, hasAttributes])

  return <div key={name}>
    {
      filtRecords.length > 0 ?
      <DropdownPieceRct
        name={name}
        changeFxn={changeFxn}
        label={label}
        value={value}
        options_in={filtRecords}
      /> :
      <DisabledTextbox
        name={name}
        label={label}
        text='Loading OR No Available Records'
      />
    }
  </div>
}
