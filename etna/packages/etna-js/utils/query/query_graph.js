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
exports.QueryGraph = void 0;
// Construct a query graph from the Magma Models,
//   so we can traverse and ask it questions.
var directed_graph_1 = require("etna-js/utils/directed_graph");
var QueryGraph = /** @class */ (function () {
    function QueryGraph(magmaModels) {
        var _this = this;
        this.unallowedModels = ['project'];
        this.includedLinkTypes = ['link', 'table'];
        this.initialized = false;
        this.models = magmaModels;
        this.graph = new directed_graph_1.DirectedGraph();
        this.allowedModels = new Set();
        // We ignore the project model and any links
        //   to project, for querying purposes only.
        // Also include linked models, to capture those paths.
        Object.entries(magmaModels).forEach(function (_a) {
            var _b = __read(_a, 2), modelName = _b[0], modelDefinition = _b[1];
            if (_this.unallowedModels.includes(modelName))
                return;
            var template = modelDefinition.template;
            _this.allowedModels.add(modelName);
            if (!_this.unallowedModels.includes(template.parent)) {
                _this.graph.addConnection(template.parent, modelName);
            }
            Object.values(template.attributes)
                .filter(function (attr) {
                return _this.includedLinkTypes.includes(attr.attribute_type) &&
                    attr.link_model_name;
            })
                .forEach(function (link) {
                if (link.link_model_name) {
                    _this.allowedModels.add(link.link_model_name);
                    _this.graph.addConnection(modelName, link.link_model_name);
                }
            });
        });
        this.initialized = true;
    }
    QueryGraph.prototype.template = function (modelName) {
        if (!Object.keys(this.models).includes(modelName))
            return null;
        return this.models[modelName].template;
    };
    QueryGraph.prototype.parentRelationship = function (modelName) {
        if (!Object.keys(this.models).includes(modelName))
            return null;
        var parentModelName = this.models[modelName].template.parent;
        if (!Object.keys(this.models).includes(parentModelName))
            return null;
        return this.models[parentModelName].template.attributes[modelName]
            .attribute_type;
    };
    QueryGraph.prototype.pathsFrom = function (modelName) {
        // this.allowedModels could have disconnected models that
        //   where only connected to models in the unallowedModels list,
        //   so they won't appear in the graph, but the user may query on
        //   them and we have to account for that.
        // An immediate example is "document".
        if (!Object.keys(this.graph.children).includes(modelName))
            return [];
        return this.graph.pathsFrom(modelName);
    };
    QueryGraph.prototype.asNormalizedHash = function (modelName) {
        // this.allowedModels could have disconnected models that
        //   where only connected to models in the unallowedModels list,
        //   so they won't appear in the graph, but the user may query on
        //   them and we have to account for that.
        // An immediate example is "document".
        if (!Object.keys(this.graph.children).includes(modelName))
            return {};
        return this.graph.asNormalizedHash(modelName);
    };
    // Here we calculate parent paths as separate entities, instead
    //   of allowing them to be in a single, flattened array.
    QueryGraph.prototype.parentPaths = function (modelName) {
        var _this = this;
        if (!Object.keys(this.graph.parents).includes(modelName))
            return [];
        var results = [];
        Object.keys(this.graph.parents[modelName]).forEach(function (p) {
            if (p !== modelName) {
                var innerPaths = _this.parentPaths(p);
                if (innerPaths.length === 0) {
                    results.push([p]);
                }
                else {
                    innerPaths.forEach(function (parentPath) {
                        parentPath.unshift(p);
                        results.push(parentPath);
                    });
                }
            }
        });
        return results;
    };
    QueryGraph.prototype.allPaths = function (modelName) {
        var _this = this;
        if (!modelName)
            return [];
        if (!Object.keys(this.graph.children).includes(modelName))
            return [];
        var parentPaths = this.parentPaths(modelName);
        // Any model that you can traverse to from any parent should
        //   also count as a path.
        // Children paths
        return this.pathsFrom(modelName)
            .concat(
        // paths up the tree
        parentPaths)
            .concat(
        // paths routing up through parents then down
        parentPaths
            .map(function (parentPath) {
            return parentPath
                .map(function (p, index) {
                return _this.pathsFrom(p).map(function (path) {
                    return parentPath.slice(0, index).concat(path);
                });
            })
                .flat(1);
        })
            .flat(1));
    };
    QueryGraph.prototype.modelHasAttribute = function (modelName, attributeName) {
        if (!this.models[modelName])
            return false;
        return !!this.models[modelName].template.attributes[attributeName];
    };
    QueryGraph.prototype.stepIsOneToMany = function (start, end) {
        // For a single model relationship (start -> end),
        //   returns `true` if it is a one-to-many
        //   relationship.
        if (!this.modelHasAttribute(start, end))
            return false;
        return ['table', 'collection'].includes(this.models[start].template.attributes[end].attribute_type);
    };
    QueryGraph.prototype.attributeIsFile = function (modelName, attributeName) {
        if (!this.modelHasAttribute(modelName, attributeName))
            return false;
        return ['file', 'image', 'file_collection'].includes(this.models[modelName].template.attributes[attributeName].attribute_type);
    };
    QueryGraph.prototype.shortestPath = function (rootModel, targetModel) {
        var paths = this.allPaths(rootModel).filter(function (potentialPath) {
            return potentialPath.includes(targetModel);
        });
        if (0 === paths.length)
            return;
        // Calculate steps to targetModel for each path
        var numberOfSteps = paths.map(function (path) {
            return path.indexOf(targetModel);
        });
        var path = paths[numberOfSteps.indexOf(Math.min.apply(Math, __spreadArray([], __read(numberOfSteps))))];
        // Direct children paths include the root, and
        //   we'll filter it out so all paths do not
        //   include the root model (eliminate redundancy).
        var pathWithoutRoot = path === null || path === void 0 ? void 0 : path.slice(0, path.indexOf(targetModel) + 1).filter(function (m) { return m !== rootModel; });
        return pathWithoutRoot;
    };
    QueryGraph.prototype.sliceableModelNamesInPath = function (startModel, endModel) {
        var _this = this;
        var modelsInPath = this.shortestPath(startModel, endModel);
        var previousModelName = startModel;
        var selectableModels = [];
        modelsInPath === null || modelsInPath === void 0 ? void 0 : modelsInPath.forEach(function (modelName) {
            if (_this.stepIsOneToMany(previousModelName, modelName)) {
                selectableModels.push(modelName);
            }
            previousModelName = modelName;
        });
        return selectableModels;
    };
    QueryGraph.prototype.childrenMap = function (modelName) {
        var _a;
        var _this = this;
        // Returns the child models, with a boolean flag
        //   for one-to-many relationships. Includes original model for convenience.
        var results = (_a = {},
            _a[modelName] = false,
            _a);
        Object.keys(this.graph.children[modelName] || {}).forEach(function (childNeighbor) {
            results[childNeighbor] = _this.stepIsOneToMany(modelName, childNeighbor);
        });
        return results;
    };
    return QueryGraph;
}());
exports.QueryGraph = QueryGraph;
