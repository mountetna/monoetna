"use strict";
var __read = (this && this.__read) || function (o, n) {
    var m = typeof Symbol === "function" && o[Symbol.iterator];
    if (!m) return o;
    var i = m.call(o), r, ar = [], e;
    try {
        while ((n === void 0 || n-- > 0) && !(r = i.next()).done) ar.push(r.value);
    }
    catch (error) { e = { error: error }; }
    finally {
        try {
            if (r && !r.done && (m = i["return"])) m.call(i);
        }
        finally { if (e) throw e.error; }
    }
    return ar;
};
var __spreadArray = (this && this.__spreadArray) || function (to, from) {
    for (var i = 0, il = from.length, j = to.length; i < il; i++, j++)
        to[j] = from[i];
    return to;
};
exports.__esModule = true;
exports.dataFrameJsonToNestedArray = exports.nestedArrayToDataFrameJson = exports.merge = exports.zipJsonDF = exports.zipDF = exports.dimensions = void 0;
var hyperformula_1 = require("hyperformula");
var maybe_1 = require("../selectors/maybe");
var dimensions = function (nestedArray) {
    return {
        numRows: nestedArray.length,
        numCols: nestedArray[0].length
    };
};
exports.dimensions = dimensions;
var zipDF = function (_a) {
    var original = _a.original, user = _a.user;
    var hfOptions = {
        licenseKey: 'gpl-v3'
    };
    // Start with the userDF, which is a subset of the original plus
    //   any added columns. Assume added columns are on the right (will have
    //   to enforce in the UI).
    // Add in missing data from the original rows.
    // And then apply equations in extra user-defined columns to
    //   ranges with hyperformula.getFillRangeData().
    var mergedDF = exports.merge({ original: original, user: user });
    var dataFrame = hyperformula_1.HyperFormula.buildFromArray(mergedDF, hfOptions);
    var userDimensions = exports.dimensions(user);
    var originalDimensions = exports.dimensions(original);
    // Extend equations and grab new cell values
    var userTopLeft = {
        sheet: 0,
        row: 1,
        col: originalDimensions.numCols
    };
    // We just grab one row to get the pattern, so we don't expect
    //   the user to populate too many of the formula cells.
    var userBottomRight = {
        sheet: 0,
        row: 1,
        col: userDimensions.numCols - 1
    };
    var newCellValues = dataFrame.getFillRangeData({
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
    return {
        formulas: dataFrame.getSheetSerialized(0),
        values: dataFrame.getSheetValues(0)
    };
};
exports.zipDF = zipDF;
var zipJsonDF = function (_a) {
    var original = _a.original, user = _a.user;
    var originalData = exports.dataFrameJsonToNestedArray(maybe_1.some(original));
    var userData = exports.dataFrameJsonToNestedArray(maybe_1.some(user));
    var values = exports.zipDF({ original: originalData, user: userData }).values;
    return exports.nestedArrayToDataFrameJson(values);
};
exports.zipJsonDF = zipJsonDF;
var merge = function (_a) {
    var original = _a.original, user = _a.user;
    var userDimensions = exports.dimensions(user);
    var originalDimensions = exports.dimensions(original);
    var numUserRows = userDimensions.numRows;
    var numUserColumns = userDimensions.numCols;
    var numOriginalColumns = originalDimensions.numCols;
    var numNewColumns = numUserColumns - numOriginalColumns;
    return original.map(function (row, index) {
        if (index < numUserRows) {
            return __spreadArray([], __read(user[index]));
        }
        else {
            // This is beyond the user preview, so
            //   we just copy the original but pad out with null values.
            return __spreadArray([], __read(row)).concat(new Array(numNewColumns).fill(null));
        }
    });
};
exports.merge = merge;
var nestedArrayToDataFrameJson = function (input) {
    var headers = input[0];
    var payload = headers.reduce(function (acc, header) {
        acc[header] = {};
        return acc;
    }, {});
    return input.slice(1).reduce(function (acc, values, rowIndex) {
        values.forEach(function (value, index) {
            var header = headers[index];
            acc[header][rowIndex.toString()] = value;
        });
        return acc;
    }, payload);
};
exports.nestedArrayToDataFrameJson = nestedArrayToDataFrameJson;
var dataFrameJsonToNestedArray = function (input) {
    if (!maybe_1.isSome(input))
        return [[]];
    var inner = Array.isArray(input) ? maybe_1.withDefault(input, {}) : input;
    var numColumns = Object.keys(inner).length;
    if (numColumns === 0)
        return [[]];
    var numRows = Object.keys(Object.values(inner)[0]).length;
    // Assume the input data is well-formed and rectangular.
    return Object.entries(inner).reduce(function (acc, _a) {
        var _b = __read(_a, 2), columnHeading = _b[0], rowData = _b[1];
        if (Object.keys(rowData).length !== numRows) {
            throw new Error('Input data is malformed and not rectangular');
        }
        acc[0].push(columnHeading);
        for (var i = 0; i < numRows; i++) {
            acc[i + 1].push(rowData[i.toString()]);
        }
        return acc;
    }, __spreadArray([], __read(new Array(1 + numRows))).map(function () { return []; }));
};
exports.dataFrameJsonToNestedArray = dataFrameJsonToNestedArray;
