"use strict";
exports.__esModule = true;
var query_selector_1 = require("../../selectors/query/query_selector");
var QuerySimplePathBuilder = /** @class */ (function () {
    function QuerySimplePathBuilder(path, rootModelName, models, flatten) {
        this.path = path;
        this.models = models;
        this.flatten = flatten;
        this.rootModelName = rootModelName;
    }
    QuerySimplePathBuilder.prototype.reducerVerb = function () {
        return this.flatten ? '::first' : '::all';
    };
    QuerySimplePathBuilder.prototype.build = function () {
        var _this = this;
        var updatedPath = [];
        var previousModelName = this.rootModelName;
        this.path.forEach(function (modelName) {
            if (query_selector_1.stepIsOneToMany(_this.models, previousModelName, modelName)) {
                updatedPath.push(modelName);
                updatedPath.push(_this.reducerVerb());
            }
            else {
                updatedPath.push(modelName);
            }
            previousModelName = modelName;
        });
        return updatedPath;
    };
    return QuerySimplePathBuilder;
}());
exports["default"] = QuerySimplePathBuilder;
