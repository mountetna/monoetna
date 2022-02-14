import React, {useMemo, useRef} from 'react';
import {HotTable} from '@handsontable/react';
import {HyperFormula} from 'hyperformula';
import Paper from '@material-ui/core/Paper';
import EditIcon from '@material-ui/icons/Edit';
import SaveIcon from '@material-ui/icons/Save';
import CancelIcon from '@material-ui/icons/Cancel';
import Button from '@material-ui/core/Button';
import {makeStyles} from '@material-ui/core/styles';

import {useModal} from 'etna-js/components/ModalDialogContainer';
import {joinNesting} from './monoids';
import {WithInputParams, DataEnvelope} from './input_types';
import {useSetsDefault} from './useSetsDefault';
import {some, Maybe} from '../../../../selectors/maybe';
import {useMemoized} from '../../../../selectors/workflow_selectors';

const useStyles = makeStyles((theme) => ({
  dialog: {
    width: '90vw',
    height: '90vh'
  }
}));

function DataTransformationModal({
  data,
  onChange
}: {
  data: any[][];
  onChange: (data: Maybe<{[key: string]: any}>) => void;
}) {
  const classes = useStyles();
  // TODO: For some reason calling up this Modal resets the Material UI
  //   theme???

  const hotTableComponent = useRef<any>(null);
  const {dismissModal} = useModal();

  const hyperformulaInstance = HyperFormula.buildEmpty({
    // to use an external HyperFormula instance,
    // initialize it with the `'internal-use-in-handsontable'` license key
    licenseKey: 'internal-use-in-handsontable'
  });

  return (
    <Paper className={classes.dialog}>
      <HotTable
        ref={hotTableComponent}
        settings={{
          data: data,
          colHeaders: true,
          rowHeaders: true,
          height: 'auto',
          licenseKey: 'non-commercial-and-evaluation',
          formulas: {
            engine: hyperformulaInstance
          },
          contextMenu: {
            items: {
              col_left: {
                name: 'Insert column to left',
                disabled: () =>
                  hotTableComponent.current.hotInstance.getSelectedLast()[1] ==
                  0
              },
              col_right: {
                name: 'Insert column to right',
                disabled: false
              },
              remove_col: {},
              undo: {},
              redo: {},
              clear_column: {}
            }
          }
        }}
      />
      <Button
        onClick={dismissModal}
        startIcon={<CancelIcon />}
        color='secondary'
        variant='contained'
      >
        Cancel
      </Button>
      <Button
        onClick={() => {
          if (hotTableComponent.current) {
            // hotTableComponent.current.hotInstance.getSourceData()
            //   returns the raw formulas, unrendered, which we could
            //   set as input values to preserve any transformations.
            const sourceData = hotTableComponent.current.hotInstance.getSourceData();
            const data = hotTableComponent.current.hotInstance.getData();

            onChange(
              some({
                source_data: toJson(sourceData),
                data: toJson(data)
              })
            );
          }

          dismissModal();
        }}
        startIcon={<SaveIcon />}
        color='primary'
        variant='contained'
      >
        Save
      </Button>
    </Paper>
  );
}

export function toJson(input: any[][]): DataEnvelope<{[key: string]: any}> {
  const headers = input[0];
  let payload = headers.reduce((acc, header) => {
    acc[header] = {};

    return acc;
  }, {});

  return input.slice(1).reduce((acc, values, rowIndex) => {
    values.forEach((value, index) => {
      let header = headers[index];
      acc[header][rowIndex.toString()] = value;
    });

    return acc;
  }, payload);
}

export function toNestedArray(
  input: DataEnvelope<{[key: string]: any}>
): any[][] {
  if (!input) return [[]];

  const numColumns = Object.keys(input).length;
  if (numColumns === 0) return [[]];

  const numRows = Object.keys(Object.values(input)[0]).length;

  // Assume the input data is well-formed and rectangular.
  return Object.entries(input).reduce(
    (
      acc: any[][],
      [columnHeading, rowData]: [string, {[key: string]: any}]
    ) => {
      if (Object.keys(rowData).length !== numRows) {
        throw new Error('Input data is malformed and not rectangular');
      }

      acc[0].push(columnHeading);

      for (var i = 0; i < numRows; i++) {
        acc[i + 1].push(rowData[i.toString()]);
      }

      return acc;
    },
    [...new Array(1 + numRows)].map(() => [])
  );
}

export default function DataTransformationInput({
  onChange,
  data,
  ...props
}: WithInputParams<
  {label?: string},
  {[key: string]: any},
  {[key: string]: any}
>) {
  const originalData = useMemoized(joinNesting, data);

  if (!originalData) return <div>No data frame!</div>;

  const {openModal} = useModal();

  const value = toNestedArray(
    useSetsDefault(
      {
        data: originalData,
        source_data: originalData
      },
      props.value,
      onChange
    ).source_data
  );

  return (
    <>
      <div>
        Your data frame has {value.length} rows and {value[0].length} columns.
      </div>
      <Button
        onClick={() => {
          openModal(
            <DataTransformationModal data={value} onChange={onChange} />
          );
        }}
        color='primary'
        startIcon={<EditIcon />}
        variant='contained'
      >
        Edit data frame
      </Button>
    </>
  );
}
