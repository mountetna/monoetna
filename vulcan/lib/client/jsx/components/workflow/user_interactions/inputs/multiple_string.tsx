import React, {useState, useMemo, useEffect, useCallback, useRef} from 'react';
import * as _ from 'lodash';

import {TYPE} from '../../../../api_types';
import {InputBackendComponent} from './input_types';
import UserInput from './user_input';
import {inputValueNonEmpty} from '../../../../selectors/workflow_selectors';
import {useBufferedInputState} from "./buffered_input_state";
import StringInput from './string';
import TextInput from 'etna-js/components/inputs/text_input';

function equalKeys(hash1: {}, hash2: {}): boolean {
  return _.isEqual(Object.keys(hash1), Object.keys(hash2));
}

const MultipleStringInput: InputBackendComponent = ({
  input,
  onChange
}) => {
  const {data} = input;

  const options: {[k: string]: string} = useMemo(() => {
    if (data) {
      return Object.values(data).reduce((a, b) => ({...a, ...b}), {});
    };

    return {};
  }, [data]);

  const [labels, setLabels] = useState<{[k: string]: string}>(options);

  const updateLabel = (newValue: string, key: string, prevLabels = labels) => {
    prevLabels[key] = newValue;
    setLabels(prevLabels);
  };

  return (
    <div>
      {Object.entries(labels).map(([k, val]) => {
        return (
          <TextInput
            header={k}
            value={val || ''}
            onChange={(newValue: string) => updateLabel(newValue, k)}
          />
        );
      })}
    </div>
  );

};

export default MultipleStringInput;
