"use strict";
var __assign = (this && this.__assign) || function () {
    __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
            s = arguments[i];
            for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p))
                t[p] = s[p];
        }
        return t;
    };
    return __assign.apply(this, arguments);
};
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
exports.createFigurePayload = exports.queryPayload = exports.userColumns = exports.isIdentifierQuery = exports.queryColumnMatrixHeadings = exports.emptyQueryClauseStamp = exports.hasMatrixSlice = exports.isMatrixSlice = exports.attributeIsFile = exports.stepIsOneToMany = exports.pathToColumn = exports.getPath = exports.selectCollectionModelNames = exports.selectMatrixModelNames = exports.selectMatrixAttributes = exports.selectAllowedModelAttributes = exports.attributeIsMatrix = exports.modelHasAttribute = void 0;
var _ = require("lodash");
var query_types_1 = require("../../contexts/query/query_types");
var modelHasAttribute = function (magmaModels, modelName, attributeName) {
    if (!magmaModels[modelName])
        return false;
    return !!magmaModels[modelName].template.attributes[attributeName];
};
exports.modelHasAttribute = modelHasAttribute;
var attributeIs = function (magmaModels, modelName, attributeName, types) {
    if (!exports.modelHasAttribute(magmaModels, modelName, attributeName))
        return false;
    return types.includes(magmaModels[modelName].template.attributes[attributeName].attribute_type);
};
var attributeIsMatrix = function (magmaModels, modelName, attributeName) {
    return attributeIs(magmaModels, modelName, attributeName, ['matrix']);
};
exports.attributeIsMatrix = attributeIsMatrix;
var selectAllowedModelAttributes = function (attributes, includeChildrenModels) {
    if (includeChildrenModels === void 0) { includeChildrenModels = false; }
    // Keep "identifier" because it's useful for ::has and ::lacks
    // Don't let folks query "up" the tree, only down it.
    var unallowedAttributeTypes = ['parent', 'link'];
    if (!includeChildrenModels) {
        unallowedAttributeTypes.push('child');
        unallowedAttributeTypes.push('collection');
        unallowedAttributeTypes.push('table');
    }
    return attributes.filter(function (attr) { return !unallowedAttributeTypes.includes(attr.attribute_type); });
};
exports.selectAllowedModelAttributes = selectAllowedModelAttributes;
var selectMatrixAttributes = function (attributes, selectedAttributes) {
    var selectedAttributeNames = selectedAttributes.map(function (attr) { return attr.attribute_name; });
    return attributes.filter(function (attr) {
        return 'matrix' === attr.attribute_type &&
            selectedAttributeNames.includes(attr.attribute_name);
    });
};
exports.selectMatrixAttributes = selectMatrixAttributes;
var selectMatrixModelNames = function (magmaModels, columns) {
    return columns
        .filter(function (column) {
        return exports.attributeIsMatrix(magmaModels, column.model_name, column.attribute_name);
    })
        .map(function (column) { return column.model_name; });
};
exports.selectMatrixModelNames = selectMatrixModelNames;
var selectCollectionModelNames = function (graph, rootModelName, selectedAttributeModelNames) {
    var sliceableModelNames = new Set();
    var fullParentage = graph.graph.fullParentage(rootModelName);
    graph.allPaths(rootModelName).forEach(function (path) {
        for (var i = 0; i < path.length - 1; i++) {
            var current = path[i];
            var next = path[i + 1];
            if ((current === rootModelName && !next) || next === rootModelName) {
                continue;
            }
            else if (i === 0 &&
                graph.stepIsOneToMany(rootModelName, current) &&
                selectedAttributeModelNames.includes(current)) {
                sliceableModelNames.add(current);
            }
            else if (graph.stepIsOneToMany(current, next) &&
                !fullParentage.includes(next) &&
                selectedAttributeModelNames.includes(next)) {
                sliceableModelNames.add(next);
            }
        }
    });
    return __spreadArray([], __read(sliceableModelNames));
};
exports.selectCollectionModelNames = selectCollectionModelNames;
var getPath = function (array, heading, currentPath) {
    if (!array)
        return [];
    if (!Array.isArray(array))
        array = [array];
    var index = array.indexOf(heading);
    if (index > -1)
        return currentPath.concat([index]);
    var innerPath = [];
    array.forEach(function (ele, index) {
        if (Array.isArray(ele)) {
            var tempPath = exports.getPath(ele, heading, currentPath.concat(index));
            if (tempPath.length > 0) {
                innerPath = tempPath;
            }
        }
    });
    return innerPath.length > 0 ? innerPath : [];
};
exports.getPath = getPath;
var pathToColumn = function (array, heading, expandMatrices) {
    var indexlessHeading = heading.split('@')[0];
    var startingIndexPlusMatrixColId = heading.split('@')[1];
    if (!startingIndexPlusMatrixColId)
        return '-1';
    var fullPath = [];
    if (expandMatrices) {
        var startingIndex = parseInt(startingIndexPlusMatrixColId.split('.')[0]);
        var sliceColId = startingIndexPlusMatrixColId.split('.')[1];
        fullPath = exports.getPath(array[startingIndex], indexlessHeading, [startingIndex]);
        if (!sliceColId)
            return fullPath.length > 0 ? fullPath[0].toString() : '-1';
        // fullPath returns the path to project#model::attribute.
        //   The tuple should be [project#model::attribute, [array, of, slice, operands]]
        var pathToSliceOperands = fullPath.slice(0, -1);
        pathToSliceOperands.push(1);
        var sliceOperands = _.at(array, pathToSliceOperands.join('.'))[0];
        if (sliceOperands == null)
            return '-1';
        var sliceIndex = sliceOperands.indexOf(sliceColId);
        // Get rid of the extra [1] used to find the slice operands
        return pathToSliceOperands.slice(0, -1).concat([sliceIndex]).join('.');
    }
    else {
        var startingIndex = parseInt(startingIndexPlusMatrixColId);
        if (!Array.isArray(array[startingIndex])) {
            if (array[startingIndex] === indexlessHeading) {
                return startingIndex.toString();
            }
            else
                return '-1';
        }
        fullPath = exports.getPath(array[startingIndex], indexlessHeading, [startingIndex]);
        return fullPath.length > 0 ? fullPath[0].toString() : '-1';
    }
};
exports.pathToColumn = pathToColumn;
var stepIsOneToMany = function (magmaModels, start, end) {
    // For a single model relationship (start -> end),
    //   returns `true` if it is a one-to-many
    //   relationship.
    return attributeIs(magmaModels, start, end, ['table', 'collection']);
};
exports.stepIsOneToMany = stepIsOneToMany;
var attributeIsFile = function (magmaModels, modelName, attributeName) {
    return attributeIs(magmaModels, modelName, attributeName, [
        'file',
        'image',
        'file_collection'
    ]);
};
exports.attributeIsFile = attributeIsFile;
var isMatrixSlice = function (slice) {
    if (!slice.clause.subclauses)
        return false;
    return slice.clause.subclauses.some(function (subclause) {
        return '::slice' === subclause.operator;
    });
};
exports.isMatrixSlice = isMatrixSlice;
var hasMatrixSlice = function (column) {
    return column.slices.some(function (slice) { return exports.isMatrixSlice(slice); });
};
exports.hasMatrixSlice = hasMatrixSlice;
var emptyQueryClauseStamp = function (modelName) {
    return __assign(__assign({}, query_types_1.EmptyQueryClause), { modelName: modelName });
};
exports.emptyQueryClauseStamp = emptyQueryClauseStamp;
var queryColumnMatrixHeadings = function (column) {
    return column.slices
        .filter(function (slice) { return exports.isMatrixSlice(slice); })
        .map(function (slice) {
        var _a;
        return (
        // Given the above isMatrixSlice filter, slice.clause.subclauses
        //   should always exist, but this makes tsc happy.
        ((_a = (slice.clause.subclauses || []).find(function (subclause) {
            return '::slice' === subclause.operator;
        })) === null || _a === void 0 ? void 0 : _a.operand).split(','));
    })
        .flat()
        .filter(function (value) { return null != value; });
};
exports.queryColumnMatrixHeadings = queryColumnMatrixHeadings;
var isIdentifierQuery = function (columns) {
    return columns.length === 1;
};
exports.isIdentifierQuery = isIdentifierQuery;
var userColumns = function (columns) {
    var columnLabels = columns.map(function (_a) {
        var display_label = _a.display_label;
        return display_label;
    });
    return exports.isIdentifierQuery(columns)
        ? [columnLabels[0], columnLabels[0]]
        : columnLabels;
};
exports.userColumns = userColumns;
var queryPayload = function (_a) {
    var query = _a.query, columns = _a.columns, expandMatrices = _a.expandMatrices;
    return {
        query: query,
        user_columns: exports.userColumns(columns),
        expand_matrices: expandMatrices,
        format: 'tsv'
    };
};
exports.queryPayload = queryPayload;
var createFigurePayload = function (_a) {
    var query = _a.query, title = _a.title, workflow = _a.workflow;
    var payload = {
        title: title,
        workflow_name: workflow.name,
        inputs: {}
    };
    Object.entries(workflow.inputQueryMap || {}).forEach(function (_a) {
        var _b = __read(_a, 2), cwlInput = _b[0], inputSource = _b[1];
        if (!query.hasOwnProperty(inputSource))
            return;
        payload.inputs[cwlInput] = query[inputSource];
        if (Array.isArray(payload.inputs[cwlInput]) ||
            (typeof payload.inputs[cwlInput] === 'object' &&
                payload.inputs[cwlInput] !== null)) {
            payload.inputs[cwlInput] = JSON.stringify(payload.inputs[cwlInput]);
        }
    });
    return payload;
};
exports.createFigurePayload = createFigurePayload;
