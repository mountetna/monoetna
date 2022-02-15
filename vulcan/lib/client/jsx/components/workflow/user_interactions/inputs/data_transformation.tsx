import React, {useRef, useState} from 'react';
import {HotTable} from '@handsontable/react';
import {HyperFormula} from 'hyperformula';
import EditIcon from '@material-ui/icons/Edit';
import SaveIcon from '@material-ui/icons/Save';
import CancelIcon from '@material-ui/icons/Cancel';
import Button from '@material-ui/core/Button';
import Dialog from '@material-ui/core/Dialog';
import DialogTitle from '@material-ui/core/DialogTitle';
import DialogContent from '@material-ui/core/DialogContent';
import DialogActions from '@material-ui/core/DialogActions';
import Typography from '@material-ui/core/Typography';
import Link from '@material-ui/core/Link';
import {makeStyles} from '@material-ui/core/styles';

import {joinNesting} from './monoids';
import {WithInputParams, DataEnvelope} from './input_types';
import {useSetsDefault} from './useSetsDefault';
import {some, Maybe} from '../../../../selectors/maybe';
import {useMemoized} from '../../../../selectors/workflow_selectors';

const useStyles = makeStyles((theme) => ({
  dialog: {
    maxWidth: '90vw',
    maxHeight: '90vh'
  },
  subtitle: {}
}));

function DataTransformationModal({
  data,
  onChange,
  onClose
}: {
  data: any[][];
  onChange: (data: Maybe<{[key: string]: any}>) => void;
  onClose: () => void;
}) {
  const classes = useStyles();

  const hotTableComponent = useRef<any>(null);
  const hyperformulaInstance = HyperFormula.buildEmpty({
    // to use an external HyperFormula instance,
    // initialize it with the `'internal-use-in-handsontable'` license key
    licenseKey: 'internal-use-in-handsontable'
  });

  return (
    <>
      <DialogTitle>Transform your data</DialogTitle>
      <DialogContent className={classes.dialog}>
        <div>
          <Typography className={classes.subtitle}>
            <Link
              target='_blank'
              rel='noreferrer'
              href='https://handsontable.github.io/hyperformula/guide/built-in-functions.html'
            >
              Supported functions
            </Link>
          </Typography>
        </div>
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
                  name: 'Insert column to left'
                },
                col_right: {
                  name: 'Insert column to right'
                },
                remove_col: {},
                undo: {},
                redo: {},
                clear_column: {}
              }
            }
          }}
        />
      </DialogContent>
      <DialogActions>
        <Button
          onClick={onClose}
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
              //   returns the raw formulas, unrendered, which we
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

            onClose();
          }}
          startIcon={<SaveIcon />}
          color='primary'
          variant='contained'
        >
          Save
        </Button>
      </DialogActions>
    </>
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
  const [open, setOpen] = useState(false);
  const originalData = useMemoized(joinNesting, data);

  function handleOnClose() {
    setOpen(false);
  }

  if (!originalData) return <div>No data frame!</div>;

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
      <Dialog open={open} onClose={handleOnClose} maxWidth='xl'>
        <DataTransformationModal
          data={value}
          onChange={onChange}
          onClose={handleOnClose}
        />
      </Dialog>
      <Button
        onClick={() => {
          setOpen(true);
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
