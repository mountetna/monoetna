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
var query_any_every_helpers_1 = require("./query_any_every_helpers");
var QueryFilterPathBuilder = /** @class */ (function () {
    function QueryFilterPathBuilder(path, rootModelName, models, anyMap) {
        this.path = path;
        this.models = models;
        this.rootModelName = rootModelName;
        this.anyMap = anyMap;
    }
    QueryFilterPathBuilder.prototype.build = function () {
        var _this = this;
        var updatedPath = [];
        var previousModelName = this.rootModelName;
        var filterAnyPath = [];
        var nestedFilterIndex = 0;
        this.path.forEach(function (modelName) {
            var foldingClause = _this.anyMap && modelName in _this.anyMap
                ? _this.anyMap[modelName]
                    ? '::any'
                    : '::every'
                : '::any';
            var newValue = [modelName, foldingClause];
            if (updatedPath.length === 0) {
                updatedPath.push.apply(updatedPath, __spreadArray([], __read(newValue)));
            }
            else {
                // here we'll nest with ::any
                filterAnyPath.push(nestedFilterIndex);
                var injected = query_any_every_helpers_1.injectValueAtPath(updatedPath, filterAnyPath, newValue);
                if (injected) {
                    // Reset this so we re-index for the new, nested array.
                    nestedFilterIndex = 0;
                }
                else {
                    // the value was not injected but rather spliced inline,
                    //   so we increment the nestedFilterIndex and pop
                    //   the last value off the filterAnyPath.
                    filterAnyPath.pop();
                }
            }
            previousModelName = modelName;
            nestedFilterIndex++;
        });
        return updatedPath;
    };
    return QueryFilterPathBuilder;
}());
exports["default"] = QueryFilterPathBuilder;
