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
var _ = require("lodash");
var FilterOperator = /** @class */ (function () {
    function FilterOperator(_a) {
        var subclause = _a.subclause, isColumnFilter = _a.isColumnFilter;
        this.subclause = subclause;
        this.isColumnFilter = isColumnFilter;
    }
    FilterOperator.allOperators = function () {
        return __spreadArray([], __read(new Set(Object.values(FilterOperator.queryOperatorsByType).reduce(function (operators, textOperatorMap) {
            operators = operators.concat(Object.values(textOperatorMap));
            return operators;
        }, []))));
    };
    FilterOperator.prototype.hasOperand = function () {
        return !(FilterOperator.terminalOperators.includes(this.subclause.operator) ||
            FilterOperator.terminalInvertOperators.includes(this.subclause.operator));
    };
    FilterOperator.prototype.hasPrepopulatedOperandOptions = function () {
        return ('string' === this.subclause.attributeType &&
            '' !== this.subclause.attributeName);
    };
    FilterOperator.prototype.attributeInputType = function () {
        switch (this.subclause.attributeType) {
            case 'string':
                return 'text';
            case 'date_time':
                return 'date';
            case 'integer':
            case 'float':
            case 'number':
                return 'number';
            case 'boolean':
                return 'boolean';
            case 'matrix':
                return 'matrix';
            case 'child':
            case 'table':
            case 'collection':
                return 'collection';
            default:
                return 'text';
        }
    };
    FilterOperator.prototype.optionsForAttribute = function () {
        return this.isColumnFilter &&
            this.attributeInputType() in FilterOperator.columnOptionsByType
            ? FilterOperator.columnOptionsByType[this.attributeInputType()]
            : this.attrOptionsWithBaseOptions();
    };
    FilterOperator.prototype.attrOptionsWithBaseOptions = function () {
        return __assign(__assign({}, FilterOperator.queryOperatorsByType.base), (FilterOperator.queryOperatorsByType[this.attributeInputType()] || {}));
    };
    FilterOperator.prototype.prettify = function () {
        return _.invert(this.optionsForAttribute())[this.subclause.operator];
    };
    FilterOperator.prototype.magmify = function (newOperator) {
        return this.optionsForAttribute()[newOperator];
    };
    FilterOperator.prototype.options = function () {
        return this.optionsForAttribute();
    };
    FilterOperator.prototype.formatOperand = function (operand) {
        return operand;
    };
    FilterOperator.queryOperatorsByType = {
        base: {
            'Is present': '::has',
            'Is missing': '::lacks'
        },
        boolean: {
            'Is true': '::true',
            'Is false': '::false',
            'Is untrue': '::untrue'
        },
        number: {
            In: '::in',
            Equals: '::=',
            'Greater than': '::>',
            'Greater than or equals': '::>=',
            'Less than': '::<',
            'Less than or equals': '::<=',
            'Not equals': '::!=',
            'Not in': '::notin'
        },
        date: {
            Equals: '::=',
            'Greater than': '::>',
            'Greater than or equals': '::>=',
            'Less than': '::<',
            'Less than or equals': '::<=',
            'Not equals': '::!='
        },
        text: {
            In: '::in',
            Equals: '::equals',
            Contains: '::matches',
            Not: '::not',
            'Not in': '::notin',
            'Greater than': '::>',
            'Greater than or equals': '::>=',
            'Less than': '::<',
            'Less than or equals': '::<='
        }
    };
    FilterOperator.columnOptionsByType = {
        matrix: {
            Slice: '::slice'
        }
    };
    FilterOperator.terminalOperators = ['::true', '::false', '::untrue'];
    FilterOperator.terminalInvertOperators = ['::has', '::lacks'];
    FilterOperator.commaSeparatedOperators = ['::in', '::slice', '::notin'];
    FilterOperator.numericTypes = ['number', 'integer', 'float'];
    return FilterOperator;
}());
exports["default"] = FilterOperator;
