import * as dataflow from 'archimedes-node/functions/dataflow/index.js';
import * as dfUtils from 'archimedes-node/functions/utils/index.js';
import HyperFormula from 'hyperformula';

function dimensions(nestedArray) {
  return {
    numRows: nestedArray.length,
    numCols: nestedArray[0].length
  };
}

function merge(original, user) {
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

const originalData = dfUtils.dataFrameJsonToNestedArray(
  await dataflow.inputJson('original_data.json')
);
const userData = dfUtils.dataFrameJsonToNestedArray(
  await dataflow.inputJson('user_data.json')
);

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

const userDimensions = dimensions(user);
const originalDimensions = dimensions(original);

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

dataflow.outputJson(
  dfUtils.nestedArrayToDataFrameJson(dataFrame.getSheetValues(0)),
  'user_data.json'
);
