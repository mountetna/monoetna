"use strict";
var __values = (this && this.__values) || function(o) {
    var s = typeof Symbol === "function" && Symbol.iterator, m = s && o[s], i = 0;
    if (m) return m.call(o);
    if (o && typeof o.length === "number") return {
        next: function () {
            if (o && i >= o.length) o = void 0;
            return { value: o && o[i++], done: !o };
        }
    };
    throw new TypeError(s ? "Object is not iterable." : "Symbol.iterator is not defined.");
};
exports.__esModule = true;
exports.applySome = exports.bindSome = exports.mapSomeAsync = exports.mapSome = exports.maybeOfNullable = exports.alternate = exports.some = exports.isSome = exports.withDefault = void 0;
function withDefault(v, d) {
    if (v == null)
        return d;
    return v[0];
}
exports.withDefault = withDefault;
function isSome(a) {
    return a != null;
}
exports.isSome = isSome;
function some(v) {
    return [v];
}
exports.some = some;
function alternate(a) {
    var e_1, _a;
    var alts = [];
    for (var _i = 1; _i < arguments.length; _i++) {
        alts[_i - 1] = arguments[_i];
    }
    if (isSome(a))
        return a;
    try {
        for (var alts_1 = __values(alts), alts_1_1 = alts_1.next(); !alts_1_1.done; alts_1_1 = alts_1.next()) {
            var alt = alts_1_1.value;
            if (isSome(alt))
                return alt;
        }
    }
    catch (e_1_1) { e_1 = { error: e_1_1 }; }
    finally {
        try {
            if (alts_1_1 && !alts_1_1.done && (_a = alts_1["return"])) _a.call(alts_1);
        }
        finally { if (e_1) throw e_1.error; }
    }
    return null;
}
exports.alternate = alternate;
function maybeOfNullable(v) {
    if (v == null)
        return null;
    return [v];
}
exports.maybeOfNullable = maybeOfNullable;
function mapSome(a, f) {
    if (isSome(a)) {
        return [f(a[0])];
    }
    return null;
}
exports.mapSome = mapSome;
function mapSomeAsync(a, f) {
    if (isSome(a)) {
        return Promise.resolve(a[0]).then(function (a) { return f(a).then(some); });
    }
    return Promise.resolve(null);
}
exports.mapSomeAsync = mapSomeAsync;
function bindSome(a, f) {
    if (isSome(a)) {
        return f(a[0]);
    }
    return null;
}
exports.bindSome = bindSome;
function applySome(f, a) {
    if (isSome(f)) {
        if (isSome(a)) {
            return [f[0](a[0])];
        }
    }
    return null;
}
exports.applySome = applySome;
