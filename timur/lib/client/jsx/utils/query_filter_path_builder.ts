import * as _ from 'lodash';
import {Model, Models} from 'etna-js/models/magma-model';
import {injectValueAtPath} from './query_any_every_helpers';
import {QueryFilterAnyMap} from '../contexts/query/query_types';

export default class QueryFilterPathBuilder {
  path: string[];
  models: Models;
  rootModelName: string;
  anyMap: QueryFilterAnyMap;

  constructor(path: string[], rootModelName: string, models: Models, anyMap: QueryFilterAnyMap) {
    this.path = path;
    this.models = models;
    this.rootModelName = rootModelName;
    this.anyMap = anyMap;
  }

  build(): any[] {
    let updatedPath: any[] = [];
    let previousModelName = this.rootModelName;
    let filterAnyPath: any[] = [];
    let nestedFilterIndex: number = 0;

    this.path.forEach((modelName: string) => {
      const foldingClause = this.anyMap && modelName in this.anyMap ? this.anyMap[modelName] ? '::any' : '::every' : '::any';

      let newValue = [
        modelName, foldingClause
      ];
      if (updatedPath.length === 0) {
        updatedPath.push(...newValue);
      } else {
        // here we'll nest with ::any
        filterAnyPath.push(nestedFilterIndex);
        let injected = injectValueAtPath(updatedPath, filterAnyPath, newValue);

        if (injected) {
          // Reset this so we re-index for the new, nested array.
          nestedFilterIndex = 0;
        } else {
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
  }
}
