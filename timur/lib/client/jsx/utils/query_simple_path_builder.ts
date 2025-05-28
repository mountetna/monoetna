import {Model, Models} from 'etna-js/models/magma-model';

export default class QuerySimplePathBuilder {
  path: string[];
  models: Models;
  flatten: boolean;
  rootModelName: string;

  constructor(
    path: string[],
    rootModelName: string,
    models: Models,
    flatten: boolean
  ) {
    this.path = path;
    this.models = models;
    this.flatten = flatten;
    this.rootModelName = rootModelName;
  }

  reducerVerb() {
    return this.flatten ? '::first' : '::all';
  }

  build(): any[] {
    let updatedPath: any[] = [];
    let previousModelName = this.rootModelName;

    this.path.forEach((modelName: string) => {
      if (this.models.model(previousModelName)?.collects(modelName)) {
        updatedPath.push(modelName);
        updatedPath.push(this.reducerVerb());
      } else {
        updatedPath.push(modelName);
      }
      previousModelName = modelName;
    });

    return updatedPath;
  }
}
