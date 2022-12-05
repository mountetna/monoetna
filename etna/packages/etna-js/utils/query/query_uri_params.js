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
exports.migrateSlices = exports.migrateSubclauses = exports.subclauseFromOldClauseFormat = exports.isOldClauseFormat = void 0;
var query_types_1 = require("../../contexts/query/query_types");
var isOldClauseFormat = function (clause) {
    return (clause.attributeName && clause.attributeType && clause.operator
    // clause.operand could be empty string as a valid entry
    );
};
exports.isOldClauseFormat = isOldClauseFormat;
var subclauseFromOldClauseFormat = function (clause) {
    return {
        attributeName: clause.attributeName || '',
        attributeType: clause.attributeType || '',
        operator: clause.operator || '',
        operand: clause.operand || ''
    };
};
exports.subclauseFromOldClauseFormat = subclauseFromOldClauseFormat;
var migrateSubclauses = function (recordFilters) {
    return recordFilters.map(function (filter) {
        return __assign(__assign({}, filter), { clauses: filter.clauses.map(function (clause) {
                var subclauses = [];
                if (exports.isOldClauseFormat(clause)) {
                    subclauses.push(exports.subclauseFromOldClauseFormat(clause));
                    subclauses = subclauses.concat(clause.subclauses || []);
                }
                else {
                    subclauses = __spreadArray([], __read((clause.subclauses || [
                        __assign({}, query_types_1.EmptyQuerySubclause)
                    ])));
                }
                return {
                    modelName: clause.modelName,
                    any: clause.any,
                    subclauses: subclauses
                };
            }) });
    });
};
exports.migrateSubclauses = migrateSubclauses;
var migrateSlices = function (columns) {
    return columns.map(function (column) {
        return __assign(__assign({}, column), { slices: column.slices.map(function (slice) {
                var clause = {
                    modelName: slice.clause.modelName,
                    any: slice.clause.any
                };
                if (exports.isOldClauseFormat(slice.clause)) {
                    clause.subclauses = [exports.subclauseFromOldClauseFormat(slice.clause)];
                }
                else {
                    clause.subclauses = __spreadArray([], __read((slice.clause.subclauses || [
                        __assign({}, query_types_1.EmptyQuerySubclause)
                    ])));
                }
                return __assign(__assign({}, slice), { clause: clause });
            }) });
    });
};
exports.migrateSlices = migrateSlices;
