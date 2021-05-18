// A port of directed_graph.rb to TypeScript

export class DirectedGraph {
  children: {[key: string]: any};
  parents: {[key: string]: any};

  constructor() {
    this.children = {};
    this.parents = {};
  }

  fullParentage(n: string): string[] {
    let result: string[] = [];

    if (!this.parents[n]) return result;

    let q: string[] = Object.keys(this.parents[n]);
    let seen: Set<string> = new Set();

    while (q.length > 0) {
      let next: string = q.splice(0, 1)[0];
      if (seen.has(next)) continue;

      seen.add(next);

      result.push(next);

      q = q.concat(Object.keys(this.parents[next]));
    }

    return [...new Set(result)];
  }

  asNormalizedHash(
    root: string,
    includeRoot: boolean = true
  ): {[key: string]: string[]} {
    let q: string[] = [root];
    let result: {[key: string]: string[]} = {};

    if (includeRoot) result[root] = [];

    let seen: Set<string> = new Set();

    while (q.length > 0) {
      let next: string = q.splice(0, 1)[0];
      if (seen.has(next)) continue;
      seen.add(next);

      let parentage: string[] = this.fullParentage(next);

      Object.keys(this.children[next]).forEach((childNode) => {
        q.push(childNode);

        if (result.hasOwnProperty(next)) result[next].push(childNode);

        parentage.forEach((grandparent) => {
          if (result.hasOwnProperty(grandparent))
            result[grandparent].push(childNode);
        });

        // Depending on the graph shape, diamonds could lead to
        //    resetting of previously calculated dependencies.
        // Here we avoid resetting existing entries in `result`
        //    and instead concatenate them if they already exist.
        if (!result.hasOwnProperty(childNode)) result[childNode] = [];
        if (result.hasOwnProperty(childNode) && result.hasOwnProperty(next))
          result[next] = result[next].concat(result[childNode]);
      });
    }

    return Object.entries(result).reduce(
      (acc: {[key: string]: string[]}, [key, value]) => {
        acc[key] = [...new Set(value)];
        return acc;
      },
      {}
    );
  }

  addConnection(parent: string, child: string) {
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
  }

  serializedPathFrom(root: string, includeRoot: boolean = true): string[] {
    let seen: Set<string> = new Set();
    let result: string[] = [];

    if (includeRoot) result.push(root);
    seen.add(root);

    let pathQ: string[][] = this.pathsFrom(root, includeRoot);
    let traversables: string[] = pathQ.flat();

    while (pathQ.length > 0) {
      let nextPath: string[] = pathQ.splice(0, 1)[0];
      if (!nextPath || nextPath.length === 0) continue;

      while (nextPath.length > 0) {
        let nextN: string = nextPath.splice(0, 1)[0];
        if (!nextN) continue;
        if (seen.has(nextN)) continue;

        if (
          Object.keys(this.parents[nextN]).some((p) => {
            return !seen.has(p) && traversables.includes(p);
          })
        ) {
          nextPath.unshift(nextN);
          pathQ.push(nextPath);
          break;
        } else {
          result.push(nextN);
          seen.add(nextN);
        }
      }
    }

    return result;
  }

  pathsFrom(root: string, includeRoot: boolean = true): string[][] {
    let result: string[][] = [];

    let parentsOfMap: {[key: string]: any} = this.descendants(
      root,
      includeRoot
    );
    let seen: Set<string> = new Set();

    Object.entries(parentsOfMap)
      .sort(([k1, p1], [k2, p2]) => {
        // if p1 is longer, sort it before p2
        if (p1.length !== p2.length) return p1.length > p2.length ? -1 : 1;
        return k1 > k2 ? 1 : -1;
      })
      .forEach(([key, parents]) => {
        if (!seen.has(key)) {
          if (Object.keys(this.children[key]).length === 0) {
            result.push(parents.concat([key]));
          } else {
            Object.keys(this.children[key])
              .sort()
              .forEach((c) => {
                result.push(parents.concat([key, c]));
              });
          }
        }

        parents.forEach((p: string) => seen.add(p));
      });

    return result;
  }

  descendants(
    parent: string,
    includeRoot: boolean = true
  ): {[key: string]: any} {
    let seen: Set<string> = new Set();
    seen.add(parent);
    let q: string[] = Object.keys(this.children[parent]);
    let parentQ: string[] = Object.keys(this.parents[parent]);

    // Because this is not an acyclical graph, the definition of descendants needs to be stronger;
    // here we believe that any path that would move through --any-- parent to this child would not be considered
    // descendant, so we first find all those parents and mark them as 'seen' so that they are not traveled.
    while (parentQ.length > 0) {
      let nextParent: string | undefined = parentQ.pop();
      if (!nextParent) break;
      if (seen.has(nextParent)) continue;
      seen.add(nextParent);
      parentQ = parentQ.concat(Object.keys(this.parents[nextParent]));
    }

    let paths: {[key: string]: any} = {};
    while (q.length > 0) {
      let child: string | undefined = q.pop();
      if (!child) break;
      if (seen.has(child)) continue;
      seen.add(child);

      if (!paths.hasOwnProperty(child))
        paths[child] = includeRoot ? [parent] : [];

      let path = paths[child];

      Object.keys(this.children[child]).forEach((childChild) => {
        q.push(childChild);
        if (!paths.hasOwnProperty(childChild))
          paths[childChild] = path.concat([child]);
      });
    }

    return paths;
  }
}
