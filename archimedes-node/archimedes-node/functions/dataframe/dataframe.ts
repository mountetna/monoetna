import * as dataflow from '../dataflow/index.js';
import HyperFormula from 'hyperformula';

type NestedArrayDataFrame = any[][];
type JsonDataFrame = {[key: string]: any};

function dimensions(nestedArray: NestedArrayDataFrame) {
  return {
    numRows: nestedArray.length,
    numCols: nestedArray[0].length
  };
}

function merge(original: NestedArrayDataFrame, user: NestedArrayDataFrame) {
  const userDimensions = dimensions(user);
  const originalDimensions = dimensions(original);
  const numUserRows = userDimensions.numRows;
  const numUserColumns = userDimensions.numCols;
  const numOriginalColumns = originalDimensions.numCols;
  const numNewColumns = numUserColumns - numOriginalColumns;
  return original.map((row, index) => {
    if (index < numUserRows - 1) {
      return [...user[index]];
    } else {
      // This is beyond the user preview, so
      //   we just copy the original but pad out with empty values.
      return [...row].concat(new Array(numNewColumns).fill(null));
    }
  });
}

export const zipDF = (original: JsonDataFrame, user: JsonDataFrame) => {
  const originalData = dataFrameJsonToNestedArray(original);
  const userData = dataFrameJsonToNestedArray(user);

  const hfOptions = {
    licenseKey: 'non-commercial-research-use'
  };

  // Start with the userDF, which is a subset of the originalData plus
  //   any added columns. Assume added columns are on the right (will have
  //   to enforce in the UI).
  // Add in missing data from the originalData rows.
  // And then apply equations in extra user-defined columns to
  //   ranges with hyperformula.getFillRangeData().
  const mergedDF = merge(originalData, userData);

  const dataFrame = HyperFormula.buildFromArray(mergedDF, hfOptions);

  const userDimensions = dimensions(userData);
  const originalDimensions = dimensions(originalData);

  // Extend equations and grab new cell values
  const userTopLeft = {
    sheet: 0,
    row: 1,
    col: originalDimensions.numCols
  };
  const userBottomRight = {
    sheet: 0,
    row: userDimensions.numRows - 1,
    col: userDimensions.numCols - 1
  };
  const newCellValues = dataFrame.getFillRangeData(
    {
      start: userTopLeft,
      end: userBottomRight
    },
    {
      start: {
        sheet: 0,
        row: 1,
        col: originalDimensions.numCols
      },
      end: {
        sheet: 0,
        row: originalDimensions.numRows - 1,
        col: userDimensions.numCols - 1
      }
    }
  );

  dataFrame.setCellContents(userTopLeft, newCellValues);

  return nestedArrayToDataFrameJson(dataFrame.getSheetValues(0));
};

// Hm...slight variations of these are also used in the Vulcan UI data_transformation.tsx
//   component. There should be a better way to avoid duplication.
export const nestedArrayToDataFrameJson = (
  input: NestedArrayDataFrame
): JsonDataFrame => {
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
};

export const dataFrameJsonToNestedArray = (
  input: JsonDataFrame
): NestedArrayDataFrame => {
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
};
