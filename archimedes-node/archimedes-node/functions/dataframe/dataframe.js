import { HyperFormula } from 'hyperformula';
const dimensions = (nestedArray) => {
    return {
        numRows: nestedArray.length,
        numCols: nestedArray[0].length
    };
};
const merge = (original, user) => {
    const userDimensions = dimensions(user);
    const originalDimensions = dimensions(original);
    const numUserRows = userDimensions.numRows;
    const numUserColumns = userDimensions.numCols;
    const numOriginalColumns = originalDimensions.numCols;
    const numNewColumns = numUserColumns - numOriginalColumns;
    return original.map((row, index) => {
        if (index < numUserRows - 1) {
            return [...user[index]];
        }
        else {
            // This is beyond the user preview, so
            //   we just copy the original but pad out with null values.
            return [...row].concat(new Array(numNewColumns).fill(null));
        }
    });
};
export const zipDF = (original, user) => {
    const originalData = dataFrameJsonToNestedArray(original);
    const userData = dataFrameJsonToNestedArray(user);
    const hfOptions = {
        licenseKey: 'gpl-v3'
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
    // We just grab two rows to get the pattern, so we don't expect
    //   the user to populate too many of the formula cells.
    const userBottomRight = {
        sheet: 0,
        row: 2,
        col: userDimensions.numCols - 1
    };
    const newCellValues = dataFrame.getFillRangeData({
        start: userTopLeft,
        end: userBottomRight
    }, {
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
    });
    dataFrame.setCellContents(userTopLeft, newCellValues);
    return nestedArrayToDataFrameJson(dataFrame.getSheetValues(0));
};
// Hm...slight variations of these are also used in the Vulcan UI data_transformation.tsx
//   component. There should be a better way to avoid duplication.
const nestedArrayToDataFrameJson = (input) => {
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
const dataFrameJsonToNestedArray = (input) => {
    if (!input)
        return [[]];
    const numColumns = Object.keys(input).length;
    if (numColumns === 0)
        return [[]];
    const numRows = Object.keys(Object.values(input)[0]).length;
    // Assume the input data is well-formed and rectangular.
    return Object.entries(input).reduce((acc, [columnHeading, rowData]) => {
        if (Object.keys(rowData).length !== numRows) {
            throw new Error('Input data is malformed and not rectangular');
        }
        acc[0].push(columnHeading);
        for (var i = 0; i < numRows; i++) {
            acc[i + 1].push(rowData[i.toString()]);
        }
        return acc;
    }, [...new Array(1 + numRows)].map(() => []));
};
