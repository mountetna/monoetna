import React, {useContext} from 'react';
import Plot from 'react-plotly.js';

import ConsignmentTable from 'etna-js/plots/components/consignment/consignment_table';
import Consignment from 'etna-js/plots/models/consignment';
import Link from 'etna-js/components/link';

import {VulcanContext} from '../../../contexts/vulcan';
import {OUTPUT_COMPONENT} from '../../../models/steps';

import {
  validDataFormat,
  uiStepInputDataRaw,
  uiStepInputDataLink
} from '../../../utils/workflow';

export default function UserOutput({step}) {
  const {pathIndex, status} = useContext(VulcanContext);

  // We need to extract the data from the input source.
  let rawInputData = uiStepInputDataRaw({step, pathIndex, status});

  let stepType = step.run.split('/')[1].replace('.cwl', '');

  stepType = Object.values(OUTPUT_COMPONENT).includes(stepType)
    ? stepType
    : 'default';

  if (!validDataFormat(rawInputData, stepType)) return null;

  const OUTPUTS = {
    default: (
      <Link link={uiStepInputDataLink({step, pathIndex, status})}>
        Download data here
      </Link>
    ),
    [OUTPUT_COMPONENT.LINK]: (
      <Link link={uiStepInputDataLink({step, pathIndex, status})}>
        Download data here
      </Link>
    ),
    [OUTPUT_COMPONENT.PLOTLY]: (
      <Plot
        data={rawInputData ? rawInputData.data : null}
        layout={rawInputData ? rawInputData.layout : null}
      ></Plot>
    ),
    [OUTPUT_COMPONENT.CONSIGNMENT]: (
      <div className='consignment-view'>
        <ConsignmentTable
          consignment={new Consignment(rawInputData)}
        ></ConsignmentTable>
      </div>
    ),
    [OUTPUT_COMPONENT.RAW]: <div className='raw-view'>{rawInputData}</div>
  };

  return OUTPUTS[stepType];
}
