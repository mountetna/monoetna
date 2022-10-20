"use strict";
// A port of directed_graph.rb to TypeScript
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
exports.DirectedGraph = void 0;
var DirectedGraph = /** @class */ (function () {
    function DirectedGraph() {
        this.children = {};
        this.parents = {};
    }
    DirectedGraph.prototype.fullParentage = function (n) {
        var result = [];
        if (!this.parents[n])
            return result;
        var q = Object.keys(this.parents[n]);
        var seen = new Set();
        while (q.length > 0) {
            var next = q.splice(0, 1)[0];
            if (seen.has(next))
                continue;
            seen.add(next);
            result.push(next);
            q = q.concat(Object.keys(this.parents[next]));
        }
        return __spreadArray([], __read(new Set(result)));
    };
    DirectedGraph.prototype.asNormalizedHash = function (root, includeRoot) {
        if (includeRoot === void 0) { includeRoot = true; }
        var q = [root];
        var result = {};
        if (includeRoot)
            result[root] = [];
        var seen = new Set();
        var _loop_1 = function () {
            var next = q.splice(0, 1)[0];
            if (seen.has(next))
                return "continue";
            seen.add(next);
            var parentage = this_1.fullParentage(next);
            Object.keys(this_1.children[next]).forEach(function (childNode) {
                q.push(childNode);
                if (result.hasOwnProperty(next))
                    result[next].push(childNode);
                parentage.forEach(function (grandparent) {
                    if (result.hasOwnProperty(grandparent))
                        result[grandparent].push(childNode);
                });
                // Depending on the graph shape, diamonds could lead to
                //    resetting of previously calculated dependencies.
                // Here we avoid resetting existing entries in `result`
                //    and instead concatenate them if they already exist.
                if (!result.hasOwnProperty(childNode))
                    result[childNode] = [];
                if (result.hasOwnProperty(childNode) && result.hasOwnProperty(next))
                    result[next] = result[next].concat(result[childNode]);
            });
        };
        var this_1 = this;
        while (q.length > 0) {
            _loop_1();
        }
        return Object.entries(result).reduce(function (acc, _a) {
            var _b = __read(_a, 2), key = _b[0], value = _b[1];
            acc[key] = __spreadArray([], __read(new Set(value)));
            return acc;
        }, {});
    };
    DirectedGraph.prototype.addConnection = function (parent, child) {
        if (!this.children.hasOwnProperty(parent)) {
            this.children[parent] = {};
        }
        if (!this.children.hasOwnProperty(child)) {
            this.children[child] = {};
        }
        this.children[parent][child] = this.children[child];
        if (!this.parents.hasOwnProperty(child)) {
            this.parents[child] = {};
        }
        if (!this.parents.hasOwnProperty(parent)) {
            this.parents[parent] = {};
        }
        this.parents[child][parent] = this.parents[parent];
    };
    DirectedGraph.prototype.serializedPathFrom = function (root, includeRoot) {
        if (includeRoot === void 0) { includeRoot = true; }
        var seen = new Set();
        var result = [];
        if (includeRoot)
            result.push(root);
        seen.add(root);
        var pathQ = this.pathsFrom(root, includeRoot);
        var traversables = pathQ.flat();
        while (pathQ.length > 0) {
            var nextPath = pathQ.splice(0, 1)[0];
            if (!nextPath || nextPath.length === 0)
                continue;
            while (nextPath.length > 0) {
                var nextN = nextPath.splice(0, 1)[0];
                if (!nextN)
                    continue;
                if (seen.has(nextN))
                    continue;
                if (Object.keys(this.parents[nextN]).some(function (p) {
                    return !seen.has(p) && traversables.includes(p);
                })) {
                    nextPath.unshift(nextN);
                    pathQ.push(nextPath);
                    break;
                }
                else {
                    result.push(nextN);
                    seen.add(nextN);
                }
            }
        }
        return result;
    };
    DirectedGraph.prototype.pathsFrom = function (root, includeRoot) {
        var _this = this;
        if (includeRoot === void 0) { includeRoot = true; }
        var result = [];
        var parentsOfMap = this.descendants(root, includeRoot);
        var seen = new Set();
        Object.entries(parentsOfMap)
            .sort(function (_a, _b) {
            var _c = __read(_a, 2), k1 = _c[0], p1 = _c[1];
            var _d = __read(_b, 2), k2 = _d[0], p2 = _d[1];
            // if p1 is longer, sort it before p2
            if (p1.length !== p2.length)
                return p1.length > p2.length ? -1 : 1;
            return k1 > k2 ? 1 : -1;
        })
            .forEach(function (_a) {
            var _b = __read(_a, 2), key = _b[0], parents = _b[1];
            if (!seen.has(key)) {
                if (Object.keys(_this.children[key]).length === 0) {
                    result.push(parents.concat([key]));
                }
                else {
                    Object.keys(_this.children[key])
                        .sort()
                        .forEach(function (c) {
                        result.push(parents.concat([key, c]));
                    });
                }
            }
            parents.forEach(function (p) { return seen.add(p); });
        });
        return result;
    };
    DirectedGraph.prototype.descendants = function (parent, includeRoot) {
        if (includeRoot === void 0) { includeRoot = true; }
        var seen = new Set();
        seen.add(parent);
        var q = Object.keys(this.children[parent]);
        var parentQ = Object.keys(this.parents[parent]);
        // Because this is not an acyclical graph, the definition of descendants needs to be stronger;
        // here we believe that any path that would move through --any-- parent to this child would not be considered
        // descendant, so we first find all those parents and mark them as 'seen' so that they are not traveled.
        while (parentQ.length > 0) {
            var nextParent = parentQ.pop();
            if (!nextParent)
                break;
            if (seen.has(nextParent))
                continue;
            seen.add(nextParent);
            parentQ = parentQ.concat(Object.keys(this.parents[nextParent]));
        }
        var paths = {};
        var _loop_2 = function () {
            var child = q.pop();
            if (!child)
                return "break";
            if (seen.has(child))
                return "continue";
            seen.add(child);
            if (!paths.hasOwnProperty(child))
                paths[child] = includeRoot ? [parent] : [];
            var path = paths[child];
            Object.keys(this_2.children[child]).forEach(function (childChild) {
                q.push(childChild);
                if (!paths.hasOwnProperty(childChild))
                    paths[childChild] = path.concat([child]);
            });
        };
        var this_2 = this;
        while (q.length > 0) {
            var state_1 = _loop_2();
            if (state_1 === "break")
                break;
        }
        return paths;
    };
    return DirectedGraph;
}());
exports.DirectedGraph = DirectedGraph;
