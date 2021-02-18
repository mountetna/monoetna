import React, {useState, useContext, useEffect} from 'react';

import ListInput from 'etna-js/components/inputs/list_input';
import SelectInput from 'etna-js/components/inputs/select_input';
import Dropdown from 'etna-js/components/inputs/dropdown';
import {
  IntegerInput,
  FloatInput
} from 'etna-js/components/inputs/numeric_input';
import Toggle from 'etna-js/components/inputs/toggle';

import {TYPE} from '../../../models/steps';
import {inputType} from '../../../selectors/workflow_selector';

export default function StepInput({step}) {
  const {workflow} = useContext(VulcanContext);

  // Select the right component depending on the
  //   workflow "inputs" type definition, or
  //   figure it out from the step itself.
  let InputComponents = step.in
    .sortBy((a, b) => a.name.localeComparison(b.name))
    .map((input) => {
      let type = inputType({workflow, input});
    });

  return <div className='step-input'>User does something here.</div>;
}
