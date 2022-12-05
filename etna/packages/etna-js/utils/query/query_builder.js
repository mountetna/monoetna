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
exports.QueryBuilder = void 0;
var query_simple_path_builder_1 = require("./query_simple_path_builder");
var query_filter_path_builder_1 = require("./query_filter_path_builder");
var query_selector_1 = require("../../selectors/query/query_selector");
var query_any_every_helpers_1 = require("./query_any_every_helpers");
var query_filter_operator_1 = require("../../components/query/query_filter_operator");
var query_uri_params_1 = require("./query_uri_params");
var QueryBuilder = /** @class */ (function () {
    function QueryBuilder(graph) {
        this.recordFilters = [];
        this.columns = [];
        this.root = '';
        this.flatten = true;
        this.orRecordFilterIndices = [];
        this.graph = graph;
    }
    QueryBuilder.prototype.addRootModel = function (modelName) {
        this.root = modelName;
    };
    QueryBuilder.prototype.addColumns = function (columns) {
        this.columns = this.columns.concat(columns);
    };
    QueryBuilder.prototype.addRecordFilters = function (recordFilters) {
        this.recordFilters = recordFilters;
    };
    QueryBuilder.prototype.setFlatten = function (flat) {
        this.flatten = flat;
    };
    QueryBuilder.prototype.setOrRecordFilterIndices = function (orRecordFilterIndices) {
        this.orRecordFilterIndices = orRecordFilterIndices;
    };
    QueryBuilder.prototype.query = function () {
        return __spreadArray(__spreadArray([
            this.root
        ], __read(this.expandedOperands(this.recordFilters))), [
            '::all',
            this.expandColumns()
        ]);
    };
    QueryBuilder.prototype.count = function () {
        return __spreadArray(__spreadArray([this.root], __read(this.expandedOperands(this.recordFilters))), ['::count']);
    };
    QueryBuilder.prototype.isNumeric = function (subclause) {
        return query_filter_operator_1["default"].numericTypes.includes(subclause.attributeType);
    };
    QueryBuilder.prototype.serializeQueryClause = function (queryClause) {
        var _this = this;
        var result = [];
        var subclauses = [];
        if (!queryClause.subclauses && query_uri_params_1.isOldClauseFormat(queryClause)) {
            subclauses.push(query_uri_params_1.subclauseFromOldClauseFormat(queryClause));
        }
        else if (queryClause.subclauses) {
            subclauses = __spreadArray([], __read(queryClause.subclauses));
        }
        subclauses.forEach(function (subclause) {
            var subclauseResult = [];
            subclauseResult.push(subclause.attributeName);
            subclauseResult.push(subclause.operator);
            if (!_this.isNumeric(subclause) &&
                query_filter_operator_1["default"].commaSeparatedOperators.includes(subclause.operator)) {
                subclauseResult.push(subclause.operand.split(','));
            }
            else if (query_filter_operator_1["default"].commaSeparatedOperators.includes(subclause.operator) &&
                _this.isNumeric(subclause)) {
                subclauseResult.push(subclause.operand.split(',').map(function (o) { return parseFloat(o); }));
            }
            else if (query_filter_operator_1["default"].terminalInvertOperators.includes(subclause.operator)) {
                // invert the model and attribute names, ignore operand
                var length = subclauseResult.length;
                var tmpOperator = subclauseResult[length - 1];
                subclauseResult[length - 1] = subclauseResult[length - 2];
                subclauseResult[length - 2] = tmpOperator;
            }
            else if (query_filter_operator_1["default"].terminalOperators.includes(subclause.operator)) {
                // ignore operand
            }
            else if (_this.isNumeric(subclause)) {
                subclauseResult.push(parseFloat(subclause.operand));
            }
            else {
                subclauseResult.push(subclause.operand);
            }
            result.push(subclauseResult);
        });
        if (1 === subclauses.length) {
            return result[0];
        }
        else {
            result.splice(0, 0, '::and');
            return result;
        }
    };
    QueryBuilder.prototype.wrapQueryClause = function (filterModelName, clause) {
        var serializedClause = this.serializeQueryClause(clause);
        if (filterModelName === clause.modelName)
            return serializedClause;
        return [
            clause.modelName,
            serializedClause,
            clause.any ? '::any' : '::every'
        ];
    };
    QueryBuilder.prototype.filterWithPath = function (filter, includeModelPath) {
        var _this = this;
        if (includeModelPath === void 0) { includeModelPath = true; }
        var result = [];
        if (filter.clauses.length > 1) {
            result = __spreadArray([
                '::and'
            ], __read(filter.clauses.map(function (clause) {
                return _this.wrapQueryClause(filter.modelName, clause);
            })));
        }
        else {
            result = this.wrapQueryClause(filter.modelName, filter.clauses[0]);
        }
        var path = this.filterPathWithModelPredicates(filter);
        if (includeModelPath && undefined != path) {
            // Inject the current [attribute, operator, operand] into
            //   the deepest array, between [model, "::any"]...
            //   to get [model, [attribute, operator, operand], "::any"]
            // At this point we know we're injecting into a tuple, so
            //   construct the valueInjectionPath that way.
            var injectionPath = query_any_every_helpers_1.nextInjectionPathItem(query_selector_1.getPath(path, filter.modelName, []));
            query_any_every_helpers_1.injectValueAtPath(path, injectionPath, result);
            result = path;
        }
        return result;
    };
    QueryBuilder.prototype.expandedOperands = function (filters) {
        var _this = this;
        var expandedFilters = [];
        var andFilters = ['::and'];
        if (this.orRecordFilterIndices.length > 0) {
            var orFilters_1 = ['::or'];
            filters.forEach(function (filter, index) {
                var expandedFilter = _this.filterWithPath(filter, _this.root !== filter.modelName);
                if (_this.orRecordFilterIndices.includes(index)) {
                    orFilters_1.push(expandedFilter);
                }
                else {
                    andFilters.push(expandedFilter);
                }
            });
            andFilters.push(orFilters_1);
            expandedFilters = [andFilters];
        }
        else if (filters.length > 1) {
            andFilters = andFilters.concat(filters.map(function (filter) {
                return _this.filterWithPath(filter, _this.root !== filter.modelName);
            }));
            expandedFilters = [andFilters];
        }
        else if (filters.length > 0) {
            // At this point, filters.length === 1...
            expandedFilters = filters.map(function (filter) {
                return _this.filterWithPath(filter, _this.root !== filter.modelName);
            });
        }
        return expandedFilters;
    };
    QueryBuilder.prototype.filterPathWithModelPredicates = function (filter) {
        var pathWithoutRoot = this.graph.shortestPath(this.root, filter.modelName);
        if (!pathWithoutRoot)
            return;
        // When constructing this path for a filter,
        //   we need to nest any collection models.
        var filterBuilder = new query_filter_path_builder_1["default"](pathWithoutRoot, this.root, this.graph.models, filter.anyMap);
        return filterBuilder.build();
    };
    QueryBuilder.prototype.slicePathWithModelPredicates = function (targetModelName) {
        var pathWithoutRoot = this.graph.shortestPath(this.root, targetModelName);
        if (!pathWithoutRoot)
            return;
        var pathBuilder = new query_simple_path_builder_1["default"](pathWithoutRoot, this.root, this.graph.models, this.flatten);
        return pathBuilder.build();
    };
    // Type should be some sort of arbitrarily nested string array,
    //   but not sure how to correctly specify all the possible permutations.
    //   [
    //     ['name'],
    //     ['species'],
    //     ['labor', 'year'],
    //     ['labor', 'completed'],
    //     ['labor', 'prize', ['name', '::equals', 'Sparta'], '::first', 'value']
    //   ]
    QueryBuilder.prototype.expandColumns = function () {
        var _this = this;
        // Convert this.attributes + this.slices into the right
        //   query format. Include the path from the root model
        //   to the attributes' model.
        if (this.columns.length === 0)
            return [''];
        var initialValues = query_selector_1.isIdentifierQuery(this.columns)
            ? this.predicateWithSlice([], this.columns[0])
            : [];
        return this.columns.slice(1).reduce(function (acc, column) {
            if (column.model_name === _this.root) {
                acc.push(_this.predicateWithSlice([], column));
            }
            else {
                var path = _this.slicePathWithModelPredicates(column.model_name);
                acc.push(_this.predicateWithSlice((path || []), column));
            }
            return acc;
        }, __spreadArray([], __read(initialValues)));
    };
    QueryBuilder.prototype.predicateWithSlice = function (path, column) {
        var _this = this;
        // If there is a slice associated with this predicate, we'll
        // inject it here, before the ::first or ::all predicate.
        var matchingSlices = column.slices || [];
        var predicate = __spreadArray([], __read(path));
        var includeAttributeName = true;
        matchingSlices.forEach(function (matchingSlice) {
            if (query_selector_1.isMatrixSlice(matchingSlice)) {
                // For matrices (i.e. ::slice), we'll construct it
                //   a little differently.
                predicate = predicate.concat(_this.serializeQueryClause(matchingSlice.clause));
                // attribute name already
                // included as part of the expanded operand
                includeAttributeName = false;
            }
            else if (_this.isTableSlice(matchingSlice)) {
                // This splicing works for tables.
                // Adds in a new array for the operand before
                //   the ::first or ::all
                var sliceModelIndex = predicate.indexOf(matchingSlice.clause.modelName);
                predicate.splice(sliceModelIndex + 1, 0, _this.serializeQueryClause(matchingSlice.clause));
            }
        });
        if (includeAttributeName) {
            predicate.push.apply(predicate, __spreadArray([], __read(this.attributeNameWithPredicate(column))));
        }
        return predicate;
    };
    QueryBuilder.prototype.isTableSlice = function (slice) {
        return !query_selector_1.isMatrixSlice(slice);
    };
    QueryBuilder.prototype.attributeNameWithPredicate = function (column) {
        // Probably only used for File / Image / FileCollection attributes?
        var predicate = [column.attribute_name];
        if (query_selector_1.attributeIsFile(this.graph.models, column.model_name, column.attribute_name)) {
            predicate.push("::" + (column.predicate || 'url'));
        }
        return predicate;
    };
    return QueryBuilder;
}());
exports.QueryBuilder = QueryBuilder;
