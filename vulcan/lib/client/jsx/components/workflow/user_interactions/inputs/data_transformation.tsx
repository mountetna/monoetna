import React, {useCallback, useEffect, useMemo, useRef, useState} from 'react';
import * as _ from 'lodash';
import {HotTable} from '@handsontable/react';
import {HyperFormula} from 'hyperformula';
import EditIcon from '@material-ui/icons/Edit';
import SaveIcon from '@material-ui/icons/Save';
import CancelIcon from '@material-ui/icons/Cancel';
import RestoreIcon from '@material-ui/icons/Restore';
import Button from '@material-ui/core/Button';
import Dialog from '@material-ui/core/Dialog';
import DialogTitle from '@material-ui/core/DialogTitle';
import DialogContent from '@material-ui/core/DialogContent';
import DialogActions from '@material-ui/core/DialogActions';
import Typography from '@material-ui/core/Typography';
import Link from '@material-ui/core/Link';
import {makeStyles} from '@material-ui/core/styles';
import CircularProgress from '@material-ui/core/CircularProgress';

import {joinNesting} from './monoids';
import {WithInputParams, DataEnvelope} from './input_types';
import {useSetsDefault} from './useSetsDefault';
import {some, Maybe, withDefault, isSome} from '../../../../selectors/maybe';
import {useMemoized} from '../../../../selectors/workflow_selectors';
import useHandsonTable from './useHandsonTable';

const useStyles = makeStyles((theme) => ({
  dialog: {
    maxWidth: '90vw',
    maxHeight: '90vh'
  },
  subtitle: {display: 'inline'},
  button: {
    margin: '1rem'
  },
  loading: {
    marginLeft: '1rem'
  }
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
  const [loading, setLoading] = useState(true);
  const classes = useStyles();

  const hotTableComponent = useRef<any>(null);
  const hyperformulaInstance = useMemo(
    () =>
      HyperFormula.buildEmpty({
        // to use an external HyperFormula instance,
        // initialize it with the `'internal-use-in-handsontable'` license key
        licenseKey: 'internal-use-in-handsontable'
      }),
    []
  );

  const columnWidths = useMemo(() => {
    let widths = new Array(data[0].length).fill(150);
    widths[0] = 200;

    return widths;
  }, [data]);

  return (
    <>
      <DialogTitle>
        Transform your data (
        <Typography className={classes.subtitle}>
          <Link
            target='_blank'
            rel='noreferrer'
            href='https://handsontable.github.io/hyperformula/guide/built-in-functions.html'
          >
            Supported functions
          </Link>
        </Typography>
        )
      </DialogTitle>
      <DialogContent className={classes.dialog}>
        <HotTable
          ref={hotTableComponent}
          settings={{
            viewportRowRenderingOffset: 10,
            viewportColumnRenderingOffset: 10,
            colWidths: columnWidths,
            autoRowSize: false,
            autoColumnSize: false,
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
                  formulaic_data: some(nestedArrayToDataFrameJson(sourceData)),
                  calculated_data: some(nestedArrayToDataFrameJson(data))
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

export function nestedArrayToDataFrameJson(
  input: any[][]
): DataEnvelope<{[key: string]: any}> {
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

export function dataFrameJsonToNestedArray(
  input: Maybe<DataEnvelope<{[key: string]: any}>>
): any[][] {
  if (!isSome(input)) return [[]];

  const inner = Array.isArray(input) ? withDefault(input, {}) : input;

  const numColumns = Object.keys(inner).length;
  if (numColumns === 0) return [[]];

  const numRows = Object.keys(Object.values(inner)[0]).length;

  // Assume the input data is well-formed and rectangular.
  return Object.entries(inner).reduce(
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
  numOutputs,
  ...props
}: WithInputParams<
  {label?: string},
  {[key: string]: any},
  {[key: string]: any}
>) {
  const [open, setOpen] = useState(false);
  const [loading, setLoading] = useState(false);
  const [isTransformed, setIsTransformed] = useState(false);

  useHandsonTable();

  const classes = useStyles();
  const originalData = useMemoized(joinNesting, data);

  function handleOnClose() {
    setOpen(false);
  }

  const destructureOnChange = useCallback(
    (v: Maybe<{[key: string]: any}>) => {
      onChange(v, true);
    },
    [onChange]
  );

  const rawData = useMemo(() => {
    return {
      calculated_data: some(originalData),
      formulaic_data: some(originalData)
    };
  }, [originalData]);

  const inputValue = useSetsDefault(rawData, props.value, destructureOnChange);

  let value: any;

  if (numOutputs === 1) {
    value = inputValue;
  } else if (inputValue.hasOwnProperty('formulaic_data')) {
    // This should work the same as the next `else` block,
    //   but will leave it here to be explicit.
    value = inputValue.formulaic_data;
  } else {
    // If the workflow author has named their CWL outputs differently, we really don't
    //   know which one is the formulaic_data. Given that the
    //   UI component uses "calculated_data" and "formulaic_data", and
    //   "formulaic_data" comes second alphabetically, we first attempt to
    //   use the second inputValue key if present. If only one key,
    //   we default to the first value.
    const inputKeys = Object.keys(inputValue).sort();
    if (inputKeys.length > 1) {
      value = inputValue[inputKeys[1]];
    } else {
      value = inputValue[inputKeys[0]];
    }
  }

  useEffect(() => {
    const transformed = !_.isEqual(
      some(nestedArrayToDataFrameJson(value)),
      rawData.formulaic_data
    );
    if (transformed !== isTransformed) setIsTransformed(transformed);
  }, [value, rawData, isTransformed]);

  const handleOnRevert = useCallback(() => {
    destructureOnChange(some(rawData));
  }, [rawData, destructureOnChange]);

  value = dataFrameJsonToNestedArray(value);

  if (!originalData || value.length === 0 || value[0].length === 0)
    return <div>No data frame!</div>;

  return (
    <>
      <div>
        Your data frame has {value.length - 1} rows and {value[0].length}{' '}
        columns. You can preview or edit the data frame now, or just click
        "Commit" to accept the raw data.
      </div>
      <Dialog
        open={open}
        onClose={handleOnClose}
        maxWidth='xl'
        TransitionProps={{
          onEntered: () => {
            setLoading(false);
          }
        }}
        disableAutoFocus={true}
        disableEnforceFocus={true}
      >
        <DataTransformationModal
          data={value}
          onChange={destructureOnChange}
          onClose={handleOnClose}
        />
      </Dialog>
      <Button
        className={classes.button}
        onClick={() => {
          setOpen(true);
          setLoading(true);
        }}
        color='primary'
        disabled={loading}
        startIcon={loading ? <CircularProgress size={24} /> : <EditIcon />}
        variant='contained'
      >
        Review or edit data frame
      </Button>
      {isTransformed ? (
        <>
          <div>** You have modified the data frame. **</div>
          <Button
            className={classes.button}
            onClick={handleOnRevert}
            color='secondary'
            startIcon={<RestoreIcon />}
            variant='contained'
          >
            Revert to raw data
          </Button>
        </>
      ) : null}
    </>
  );
}
