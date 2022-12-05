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
exports.injectValueAtPath = exports.nextInjectionPathItem = void 0;
var _ = require("lodash");
var anyOrEvery = ['::any', '::every'];
var nextInjectionPathItem = function (injectionPath) {
    var nextItemPath = __spreadArray([], __read(injectionPath));
    nextItemPath[injectionPath.length - 1] =
        injectionPath[injectionPath.length - 1] + 1;
    return nextItemPath;
};
exports.nextInjectionPathItem = nextInjectionPathItem;
var injectValueAtPath = function (array, valueInjectionPath, value) {
    var currentValue = _.get(array, valueInjectionPath);
    if (anyOrEvery.includes(currentValue)) {
        _.set(array, valueInjectionPath, value);
        // We create a new injection path to inject a final "::any" or "::every"
        //   after the value.
        var anyInjectionPath = exports.nextInjectionPathItem(valueInjectionPath);
        _.set(array, anyInjectionPath, currentValue);
        return true;
    }
    else {
        // We need to "splice" in the values at the path...
        var refArray_1 = array;
        valueInjectionPath.forEach(function (index) {
            if (Array.isArray(refArray_1[index])) {
                refArray_1 = refArray_1[index];
            }
            else {
                // We are at the injection spot
                refArray_1.splice.apply(refArray_1, __spreadArray([index, 0], __read(value)));
            }
        });
        return false;
    }
};
exports.injectValueAtPath = injectValueAtPath;
