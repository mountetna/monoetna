// Construct a query graph from the Magma Models,
//   so we can traverse and ask it questions.
import {DirectedGraph} from 'etna-js/utils/directed_graph';
import {Model, Attribute, Models, ModelsObject} from 'etna-js/models/magma-model';

export class QueryGraph {
  models: Models;
  graph: DirectedGraph;
  includedLinkTypes: string[] = ['link'];
  initialized: boolean = false;

  constructor(magmaModels: ModelsObject) {
    this.models = new Models(magmaModels);
    this.graph = new DirectedGraph();

    // We ignore the project model and any links
    //   to project, for querying purposes only.
    // Also include linked models, to capture those paths.
    this.models.each(
      (modelName: string, model: Model) => {
        if (model.parent) this.graph.addConnection(model.parent, modelName);

        model.selectAttributes(
            (attr: Attribute) =>
              !!(this.includedLinkTypes.includes(attr.attribute_type) && attr.link_model_name)
          )
          .forEach((link: Attribute) => {
            this.graph.addConnection(link.link_model_name as string, modelName);
          });
      }
    );

    this.initialized = true;
  }

  connectedModelsAnd(modelName: string | null): Set<string> {
    return modelName ? this.connectedModels(modelName).add(modelName) : new Set();
  }

  connectedModels(modelName: string | null): Set<string> {
    return modelName ? new Set(this.allPaths(modelName).flat()) : new Set();
  }

  ancestors(modelName: string | null): Set<string> {
    return modelName ? new Set(this.graph.ancestors(modelName)) : new Set();
  }

  sliceable(modelName: string, selectedModel: string): boolean {
    if (!this.models.has(modelName)) return false;

    return this.allPaths(modelName).some((path: string[]) => {
      for (let i = 0; i < path.length - 1; i++) {
        let current = path[i];
        let next = path[i + 1];

        if ((current === modelName && !next) ||
          next === modelName) continue;

        if (
          i === 0 &&
          this.models.model(modelName)!.collects(current) &&
          selectedModel == current
        ) {
          return true;
        } else if (
          this.models.model(current)?.collects(next) &&
          !this.graph.ancestors(modelName).includes(next) &&
          selectedModel == next
        ) {
          return true;
        }
      }
      return false;
    });
  }

  parentRelationship(modelName: string): string | undefined {
    const model = this.models.model(modelName);

    if (!model) return undefined;

    const parent: Attribute =  model.selectAttributes(
      (a: Attribute) => a.attribute_type == 'parent'
    )[0];

    if (!parent) return undefined;

    return parent.link_attribute_type;
  }

  pathsFrom(modelName: string): string[][] {
    //   so they won't appear in the graph, but the user may query on
    //   them and we have to account for that.
    // An immediate example is "document".
    if (!Object.keys(this.graph.children).includes(modelName)) return [];
    return this.graph.pathsFrom(modelName);
  }

  asNormalizedHash(modelName: string): {[key: string]: string[]} {
    if (!Object.keys(this.graph.children).includes(modelName)) return {};
    return this.graph.asNormalizedHash(modelName);
  }

  // Here we calculate parent paths as separate entities, instead
  //   of allowing them to be in a single, flattened array.
  parentPaths(modelName: string): string[][] {
    if (!(modelName in this.graph.parents)) return [];

    let results: string[][] = [];

    Object.keys(this.graph.parents[modelName]).forEach((p: string) => {
      if (p !== modelName) {
        let innerPaths = this.parentPaths(p);
        if (innerPaths.length === 0) {
          results.push([p]);
        } else {
          innerPaths.forEach((parentPath: string[]) => {
            parentPath.unshift(p);
            results.push(parentPath);
          });
        }
      }
    });

    return results;
  }

  allPaths(modelName: string | null): string[][] {
    if (!modelName) return [];

    if (!Object.keys(this.graph.children).includes(modelName)) return [];

    let parentPaths = this.parentPaths(modelName);

    // Any model that you can traverse to from any parent should
    //   also count as a path.

    // Children paths
    return this.pathsFrom(modelName)
      .concat(
        // paths up the tree
        parentPaths
      )
      .concat(
        // paths routing up through parents then down
        parentPaths.map((parentPath: string[]) =>
          parentPath.map((p: string, index: number) =>
            p == 'project' ? [] : this.pathsFrom(p).map(
              path => path.includes(modelName) ? [] : parentPath.slice(0, index).concat(path)
            )
          ).flat(1)
        ).flat(1).filter(path => path.length > 0)
      );
  }

  shortestPath(rootModel: string, targetModel: string): string[] | undefined {
    let paths = this.allPaths(rootModel).filter((potentialPath: string[]) =>
      potentialPath.includes(targetModel)
    );

    if (0 === paths.length) return;

    // Calculate steps to targetModel for each path
    let numberOfSteps = paths.map((path: string[]) => {
      return path.indexOf(targetModel);
    });
    let path = paths[numberOfSteps.indexOf(Math.min(...numberOfSteps))];

    // Direct children paths include the root, and
    //   we'll filter it out so all paths do not
    //   include the root model (eliminate redundancy).
    const pathWithoutRoot = path
      ?.slice(0, path.indexOf(targetModel) + 1)
      .filter((m) => m !== rootModel);

    return pathWithoutRoot;
  }

  sliceableModelNamesInPath(startModel: string, endModel: string): string[] {
    let modelsInPath = this.shortestPath(startModel, endModel);
    let previousModelName = startModel;
    let selectableModels: string[] = [];

    modelsInPath?.forEach((modelName) => {
      if (this.models.model(previousModelName)?.collects(modelName)) {
        selectableModels.push(modelName);
      }
      previousModelName = modelName;
    });

    return selectableModels;
  }

  subgraphMap(modelName: string): {[key: string]: boolean} {
    return {
      ...this.childrenMap( modelName ),
      [ modelName ]: false
    };
  }

  childrenMap(modelName: string): {[key: string]: boolean} {
    // Returns the child models, with a boolean flag
    //   for one-to-many relationships. Includes original model for convenience.
    let results: {[key: string]: boolean} = { };

    Object.keys(this.graph.children[modelName] || {}).forEach(
      (childNeighbor: string) => {
        results[childNeighbor] = !!this.models.model(modelName)?.collects(childNeighbor);
      }
    );

    return results;
  }
}
