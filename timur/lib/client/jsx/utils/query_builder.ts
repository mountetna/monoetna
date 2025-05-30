import * as _ from 'lodash';
import {
  QueryColumn,
  QueryFilter,
  QuerySlice,
  QueryClause,
  QuerySubclause
} from '../contexts/query/query_types';
import {QueryGraph} from './query_graph';
import QuerySimplePathBuilder from './query_simple_path_builder';
import QueryFilterPathBuilder from './query_filter_path_builder';
import {
  isMatrixSlice,
  getPath,
  isIdentifierQuery
} from '../selectors/query_selector';
import {
  injectValueAtPath,
  nextInjectionPathItem
} from './query_any_every_helpers';
import FilterOperator from '../components/query/query_filter_operator';
import {
  isOldClauseFormat,
  subclauseFromOldClauseFormat
} from './query_uri_params';

export class QueryBuilder {
  graph: QueryGraph;
  recordFilters: QueryFilter[] = [];
  columns: QueryColumn[] = [];
  root: string = '';
  flatten: boolean = true;
  orRecordFilterIndices: number[] = [];

  constructor(graph: QueryGraph) {
    this.graph = graph;
  }

  addRootModel(modelName: string) {
    this.root = modelName;
  }

  addColumns(columns: QueryColumn[]) {
    this.columns = this.columns.concat(columns);
  }

  addRecordFilters(recordFilters: QueryFilter[]) {
    this.recordFilters = recordFilters;
  }

  setFlatten(flat: boolean) {
    this.flatten = flat;
  }

  setOrRecordFilterIndices(orRecordFilterIndices: number[]) {
    this.orRecordFilterIndices = orRecordFilterIndices;
  }

  query(): any[] {
    return [
      this.root,
      ...this.expandedOperands(this.recordFilters),
      '::all',
      this.expandColumns()
    ];
  }

  count(): any[] {
    return [this.root, ...this.expandedOperands(this.recordFilters), '::count'];
  }

  isNumeric(subclause: QuerySubclause): boolean {
    return FilterOperator.numericTypes.includes(subclause.attributeType);
  }

  serializeQueryClause(queryClause: QueryClause): any[] {
    let result: any[] = [];
    let subclauses: QuerySubclause[] = [];

    if (!queryClause.subclauses && isOldClauseFormat(queryClause)) {
      subclauses.push(subclauseFromOldClauseFormat(queryClause));
    } else if (queryClause.subclauses) {
      subclauses = [...queryClause.subclauses];
    }

    subclauses.forEach((subclause: QuerySubclause) => {
      let subclauseResult: any[] = [];

      subclauseResult.push(subclause.attributeName);
      subclauseResult.push(subclause.operator);

      if (
        !this.isNumeric(subclause) &&
        FilterOperator.commaSeparatedOperators.includes(subclause.operator)
      ) {
        subclauseResult.push((subclause.operand as string).split(','));
      } else if (
        FilterOperator.commaSeparatedOperators.includes(subclause.operator) &&
        this.isNumeric(subclause)
      ) {
        subclauseResult.push(
          (subclause.operand as string).split(',').map((o) => parseFloat(o))
        );
      } else if (
        FilterOperator.terminalInvertOperators.includes(subclause.operator)
      ) {
        // invert the model and attribute names, ignore operand
        let length = subclauseResult.length;
        let tmpOperator = subclauseResult[length - 1];
        subclauseResult[length - 1] = subclauseResult[length - 2];
        subclauseResult[length - 2] = tmpOperator;
      } else if (
        FilterOperator.terminalOperators.includes(subclause.operator)
      ) {
        // ignore operand
      } else if (this.isNumeric(subclause)) {
        subclauseResult.push(parseFloat(subclause.operand as string));
      } else {
        subclauseResult.push(subclause.operand);
      }

      result.push(subclauseResult);
    });

    if (1 === subclauses.length) {
      return result[0];
    } else {
      result.splice(0, 0, '::and');
      return result;
    }
  }

  wrapQueryClause(filterModelName: string, clause: QueryClause): any[] {
    const serializedClause = this.serializeQueryClause(clause);

    if (filterModelName === clause.modelName) return serializedClause;

    return [
      clause.modelName,
      serializedClause,
      clause.any ? '::any' : '::every'
    ];
  }

  filterWithPath(filter: QueryFilter, includeModelPath: boolean = true): any[] {
    let result: any[] = [];

    if (filter.clauses.length > 1) {
      result = [
        '::and',
        ...filter.clauses.map((clause) =>
          this.wrapQueryClause(filter.modelName, clause)
        )
      ];
    } else {
      result = this.wrapQueryClause(filter.modelName, filter.clauses[0]);
    }

    let path: string[] | undefined = this.filterPathWithModelPredicates(filter);
    if (includeModelPath && undefined != path) {
      // Inject the current [attribute, operator, operand] into
      //   the deepest array, between [model, "::any"]...
      //   to get [model, [attribute, operator, operand], "::any"]
      // At this point we know we're injecting into a tuple, so
      //   construct the valueInjectionPath that way.
      let injectionPath = nextInjectionPathItem(
        getPath(path, filter.modelName, [])
      );
      injectValueAtPath(path, injectionPath, result);
      result = path;
    }

    return result;
  }

  expandedOperands(filters: QueryFilter[]) {
    let expandedFilters: any[] = [];
    let andFilters: any[] = ['::and'];

    if (this.orRecordFilterIndices.length > 0) {
      let orFilters: any[] = ['::or'];

      filters.forEach((filter, index: number) => {
        let expandedFilter = this.filterWithPath(
          filter,
          this.root !== filter.modelName
        );

        if (this.orRecordFilterIndices.includes(index)) {
          orFilters.push(expandedFilter);
        } else {
          andFilters.push(expandedFilter);
        }
      });

      andFilters.push(orFilters);
      expandedFilters = [andFilters];
    } else if (filters.length > 1) {
      andFilters = andFilters.concat(
        filters.map((filter) =>
          this.filterWithPath(filter, this.root !== filter.modelName)
        )
      );
      expandedFilters = [andFilters];
    } else if (filters.length > 0) {
      // At this point, filters.length === 1...
      expandedFilters = filters.map((filter) =>
        this.filterWithPath(filter, this.root !== filter.modelName)
      );
    }
    return expandedFilters;
  }

  filterPathWithModelPredicates(filter: QueryFilter): any[] | undefined {
    const pathWithoutRoot = this.graph.shortestPath(
      this.root,
      filter.modelName
    );

    if (!pathWithoutRoot) return;

    // When constructing this path for a filter,
    //   we need to nest any collection models.
    const filterBuilder = new QueryFilterPathBuilder(
      pathWithoutRoot,
      this.root,
      this.graph.models,
      filter.anyMap
    );
    return filterBuilder.build();
  }

  slicePathWithModelPredicates(targetModelName: string): any[] | undefined {
    const pathWithoutRoot = this.graph.shortestPath(this.root, targetModelName);

    if (!pathWithoutRoot) return;

    const pathBuilder = new QuerySimplePathBuilder(
      pathWithoutRoot,
      this.root,
      this.graph.models,
      this.flatten
    );
    return pathBuilder.build();
  }

  // Type should be some sort of arbitrarily nested string array,
  //   but not sure how to correctly specify all the possible permutations.
  //   [
  //     ['name'],
  //     ['species'],
  //     ['labor', 'year'],
  //     ['labor', 'completed'],
  //     ['labor', 'prize', ['name', '::equals', 'Sparta'], '::first', 'value']
  //   ]
  expandColumns(): (string | string[] | (string | string[])[])[] {
    // Convert this.attributes + this.slices into the right
    //   query format. Include the path from the root model
    //   to the attributes' model.
    if (this.columns.length === 0) return [''];

    let initialValues = isIdentifierQuery(this.columns)
      ? this.predicateWithSlice([], this.columns[0])
      : [];

    return this.columns.slice(1).reduce(
      (acc: any[], column: QueryColumn) => {
        if (column.model_name === this.root) {
          acc.push(this.predicateWithSlice([], column));
        } else {
          let path = this.slicePathWithModelPredicates(column.model_name);

          acc.push(
            this.predicateWithSlice((path || []) as string[], column) as (
              | string
              | string[]
            )[]
          );
        }

        return acc;
      },
      [...initialValues]
    );
  }

  predicateWithSlice(
    path: string[],
    column: QueryColumn
  ): (string | string[] | (string | string[] | number)[])[] {
    // If there is a slice associated with this predicate, we'll
    // inject it here, before the ::first or ::all predicate.
    let matchingSlices = column.slices || [];

    let predicate: (string | string[] | (string | string[] | number)[])[] = [
      ...path
    ];

    let includeAttributeName = true;

    matchingSlices.forEach((matchingSlice: QuerySlice) => {
      if (isMatrixSlice(matchingSlice)) {
        // For matrices (i.e. ::slice), we'll construct it
        //   a little differently.
        predicate = predicate.concat(
          this.serializeQueryClause(matchingSlice.clause)
        );
        // attribute name already
        // included as part of the expanded operand
        includeAttributeName = false;
      } else if (this.isTableSlice(matchingSlice)) {
        // This splicing works for tables.
        // Adds in a new array for the operand before
        //   the ::first or ::all
        let sliceModelIndex = predicate.indexOf(matchingSlice.clause.modelName);
        predicate.splice(
          sliceModelIndex + 1,
          0,
          this.serializeQueryClause(matchingSlice.clause)
        );
      }
    });

    if (includeAttributeName)
      {predicate.push(...this.attributeNameWithPredicate(column));}

    return predicate;
  }

  isTableSlice(slice: QuerySlice) {
    return !isMatrixSlice(slice);
  }

  attributeNameWithPredicate(column: QueryColumn) {
    // Probably only used for File / Image / FileCollection attributes?
    let predicate = [column.attribute_name];
    if (
      this.graph.models.attribute(
	column.model_name, column.attribute_name
      )?.isFile()
    ) {
      predicate.push(`::${column.predicate || 'url'}`);
    }

    return predicate;
  }
}
