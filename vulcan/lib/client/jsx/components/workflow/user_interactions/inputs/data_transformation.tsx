import React, {useCallback, useEffect, useMemo, useRef, useState} from 'react';
import * as _ from 'lodash';
import {HotTable} from '@handsontable/react';
import {HyperFormula} from 'hyperformula';
import EditIcon from '@material-ui/icons/Edit';
import TableChartIcon from '@material-ui/icons/TableChart';
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
import {WithInputParams} from './input_types';
import {useSetsDefault} from './useSetsDefault';
import {some, Maybe} from '../../../../selectors/maybe';
import {useMemoized} from '../../../../selectors/workflow_selectors';
import useHandsonTable from './useHandsonTable';

import {
  zipDF,
  dimensions,
  merge,
  dataFrameJsonToNestedArray,
  nestedArrayToDataFrameJson,
  NestedArrayDataFrame,
  JsonDataFrame
} from 'etna-js/utils/dataframe';

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
  },
  helpdoc: {
    maxWidth: '600px',
    marginTop: '1rem',
    marginBottom: '1rem'
  },
  propagateButton: {
    marginBottom: '1rem'
  }
}));

function DataTransformationModal({
  userData,
  originalData,
  onChange,
  onClose,
  dialogTitle,
  dialogTexts
}: {
  userData: NestedArrayDataFrame;
  originalData: NestedArrayDataFrame;
  onChange: (data: Maybe<JsonDataFrame>) => void;
  onClose: () => void;
  dialogTitle: string;
  dialogTexts: string[];
}) {
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

  const numOriginalCols = dimensions(originalData).numCols;

  const columnWidths = useMemo(() => {
    let widths = new Array(dimensions(userData).numCols).fill(150);
    widths[0] = 200;

    return widths;
  }, [userData]);

  const isReadOnlyColumn = useCallback(
    (colIndex: number) => {
      return colIndex < numOriginalCols;
    },
    [numOriginalCols]
  );

  const canInsertColumn = useCallback(
    (colIndex: number) => {
      return colIndex < numOriginalCols - 1;
    },
    [numOriginalCols]
  );

  const handleExtendFormulas = useCallback(() => {
    if (!hotTableComponent.current) return;

    const zippedData = zipDF({
      original: originalData,
      user: hotTableComponent.current.hotInstance.getSourceData()
    });

    hotTableComponent.current.hotInstance.loadData(zippedData.formulas);
  }, [originalData, hotTableComponent]);

  const mergedData = useMemo(() => {
    return merge({original: originalData, user: userData});
  }, [userData, originalData]);

  return (
    <>
      <DialogTitle>
        {dialogTitle} (
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
        {
          dialogTexts==null ? null : dialogTexts.map(
            (val) => <Typography className={classes.helpdoc}>{val}</Typography>
          )
        }
        <Button
          className={classes.propagateButton}
          onClick={handleExtendFormulas}
          startIcon={<TableChartIcon />}
          color='secondary'
          variant='contained'
        >
          Propagate Formulas
        </Button>
        <HotTable
          ref={hotTableComponent}
          settings={{
            viewportRowRenderingOffset: 10,
            viewportColumnRenderingOffset: 10,
            colWidths: columnWidths,
            autoRowSize: false,
            autoColumnSize: false,
            data: mergedData,
            colHeaders: true,
            rowHeaders: true,
            height: 'auto',
            licenseKey: 'non-commercial-and-evaluation',
            formulas: {
              engine: hyperformulaInstance
            },
            cells: (row: number, col: number, props: any) => {
              if (row > 0 && col < numOriginalCols) {
                return {
                  readOnly: true
                };
              }

              return {};
            },
            contextMenu: {
              items: {
                col_right: {
                  name: 'Insert column to right',
                  disabled: () => {
                    return canInsertColumn(
                      hotTableComponent.current.hotInstance.getSelectedLast()[1]
                    );
                  }
                },
                remove_col: {
                  disabled: () => {
                    return isReadOnlyColumn(
                      hotTableComponent.current.hotInstance.getSelectedLast()[1]
                    );
                  }
                },
                undo: {},
                redo: {},
                clear_column: {
                  disabled: () => {
                    return isReadOnlyColumn(
                      hotTableComponent.current.hotInstance.getSelectedLast()[1]
                    );
                  }
                }
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

              // Only send back the first couple rows of data,
              //   since the server will extend formulas based on the
              //   first couple of rows only. This will make sure we
              //   aren't saving or sending giant blobs of data
              //   as an input.
              const truncatedFormulas = sourceData.slice(0, 21);
              const truncatedData = data.slice(0, 21);

              onChange(
                some({
                  formulaic_data: some(
                    nestedArrayToDataFrameJson(truncatedFormulas)
                  ),
                  calculated_data: some(
                    nestedArrayToDataFrameJson(truncatedData)
                  )
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

export default function HandsOnTableInput({
  onChange,
  data,
  numOutputs,
  ...props
}: WithInputParams<
  {label?: string},
  {[key: string]: any},
  {[key: string]: any}
>,
  dataName: string,
  afterSummaryText: string,
  openButtonText: string,
  ifChangedText: string,
  resetButtonText: string,
  dialogTitle: string,
  dialogTexts: string[],
) {
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

  const originalAsNestedArray = useMemo(() => {
    return dataFrameJsonToNestedArray(some(originalData));
  }, [originalData]);

  if (!originalData || value.length === 0 || value[0].length === 0)
    {return <div>No data frame!</div>;}

  return (
    <>
      <div>
        Your {dataName} has {originalAsNestedArray.length - 1} rows and{' '}
        {value[0].length} columns. {afterSummaryText}
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
          userData={value}
          originalData={originalAsNestedArray}
          onChange={destructureOnChange}
          onClose={handleOnClose}
          dialogTitle={dialogTitle}
          dialogTexts={dialogTexts}
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
        {openButtonText}
      </Button>
      {isTransformed ? (
        <>
          <div>{ifChangedText}</div>
          <Button
            className={classes.button}
            onClick={handleOnRevert}
            color='secondary'
            startIcon={<RestoreIcon />}
            variant='contained'
          >
            {resetButtonText}
          </Button>
        </>
      ) : null}
    </>
  );
}

export function AnnotationEditorInput({
  onChange,
  data,
  numOutputs,
  ...props
}: WithInputParams<
  {label?: string},
  {[key: string]: any},
  {[key: string]: any}
>) {
  return HandsOnTableInput({
    onChange,
    data,
    numOutputs,
    ...props
  },
  'Annotation data frame',
  'Click the button below to add annotations.',
  'Edit Annotations',
  'Data frame changed from blank original. Click below to reset',
  'Revert to original',
  'Add/Edit Annotations',
  [
    'This is your annotations data frame. Add columns to the right to record annotations or notes and use column names (top row) to establish what columns represent.',
    '- Annotations: columns with names starting with \'annot\', e.g, \'annots_fine\' or \'annotations_broad\', will be treated as cluster calls that should be pulled into the dataset.',
    '- Notes: columns whose names do not start with \'annot\' are "free" columns that will appear in your csv and xlsx downloads, but will be otherwise ignored by this workflow.',
    'To apply a formula in a column, establish the formula in the first cell (2nd row) by starting with an \'=\'. As an example, \'=IF(ISNUMBER(OR(FIND("CD4",B2),FIND(\"CD8\",B2))), \"T\", \"non-T\")\', might be useful for initiating an \'annots_broad\' column from an \'annots_fine\' column-B. Then click the "Propagate Formulas" button to propagate the formula down the rest of the column.',
    'Save and Commit often to not lose work!'
  ])
}

export function DataTransformationInput({
  onChange,
  data,
  numOutputs,
  ...props
}: WithInputParams<
  {label?: string},
  {[key: string]: any},
  {[key: string]: any}
>) {
  return HandsOnTableInput({
    onChange,
    data,
    numOutputs,
    ...props
  },
  'data frame',
  'You can preview or edit the data frame now, or just click "Commit" to accept the raw data.',
  'Review or edit data frame',
  '** You have modified the data frame. **',
  'Revert to raw data',
  'Transform your data',
  [
    'This is a preview of your data frame. You can edit the column headings or append additional columns on the right, by right-clicking and selecting "Insert column to right" in the context menu.',
    'To apply a formula in a new column, just add a couple of cells with the formula to establish the pattern. Click the "Propagate Formulas" button to propagate the formulas to the entire table. Save, Commit, and Run!'
  ])
}
