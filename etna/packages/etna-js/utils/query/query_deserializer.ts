/* Take a query string + user_columns and convert them into the JSON
  objects that we use in the Contexts.
May also be useful for the URL params, at some point, if we want to shorten them?
i.e.

query: ["patient",["sample", "tissue_type", "::=", "Primary"], "::all",["ipi_number"]]
user_columns: ["ipi_number","ipi_number"]

should become something like:

QueryColumnContext
{
     columns: [{
      model_name: "patient",
      attribute_name: "ipi_number",
      display_label: "ipi_number",  <-- note that this may not quite work for matrices...
      slices: []
    }]
}

QueryWhereContext
{
  recordFilters: [{
    model_name: "sample",
    anyMap: {
        "sample": true
    },
    clauses: [{
        any: true,
        model_name: 'sample',
        subclauses: [{
            attributeName: 'tumor_type',
            attributeType: 'string',
            operator: '::=',
            operand: 'Primary'
        }]
    }]
  }],
  orRecordFilterIndices: []
}

QueryGraphText
{
  rootModel: "patient"
}
*/
import {
  QueryFilter,
  QueryColumn,
  QueryClause,
  QuerySubclause,
  EmptyQuerySubclause
} from '../../contexts/query/query_types';
import QueryFilterOperator from '../../components/query/query_filter_operator';
export class QueryDeserializer {
  query: any[];
  userColumns: string[];
  operatorHelper: QueryFilterOperator;

  and: string;
  or: string;
  all: string;
  every: string;
  any: string;

  constructor(query: any[], userColumns: string[]) {
    this.query = query;
    this.userColumns = userColumns;

    this.and = '::and';
    this.or = '::or';
    this.all = '::all';
    this.every = '::every';
    this.any = '::any';

    this.operatorHelper = new QueryFilterOperator({
      subclause: {...EmptyQuerySubclause},
      isColumnFilter: false
    });
  }

  rootModel() {
    return this.query[0];
  }

  flattenFilters(rawFilters: any[]): QueryFilter[] {
    let filterArrays: any[] = [];
    let queue: any[] = [...rawFilters];

    while (queue.length > 0) {
      let filter = queue.shift();

      let firstArg = filter[0];
      if (filter === this.and || filter === this.or) {
        continue;
      } else if (
        Array.isArray(filter) &&
        (firstArg === this.and || firstArg === this.or)
      ) {
        filter.slice(1).forEach((nestedFilter) => {
          queue.push(nestedFilter);
        });
      } else {
        filterArrays.push(this.toFilter(filter));
      }
    }

    return filterArrays;
  }

  recordFilters(): QueryFilter[] {
    // If the first element of this.rawFilters() is "::and", then there
    //   must be more than 1 filter.
    // Otherwise, just one filter, that we return.
    const rawFilters = this.rawFilters();

    if (this.isCompositeFilter(rawFilters)) {
      // This could still be an ::and on the rootModel(), or
      //   could include paths to other models...deal with that
      //   separately.
      return this.flattenFilters(rawFilters);
    } else {
      return [this.toFilter(rawFilters)];
    }
  }

  isAndFilter(rawFilter: any[]): boolean {
    return this.and === rawFilter[0];
  }

  isOrFilter(rawFilter: any[]): boolean {
    return this.or === rawFilter[0];
  }

  isCompositeFilter(rawFilter: any[]): boolean {
    // If the first element is ::and or ::or, then is a composition of
    //   multiple filters.
    return this.isAndFilter(rawFilter) || this.isOrFilter(rawFilter);
  }

  isCollectionModelFilter(rawFilterArray: any[]): boolean {
    return [this.any, this.every].includes(
      rawFilterArray[rawFilterArray.length - 1]
    );
  }

  toFilter(rawFilterArray: any[]): QueryFilter {
    // Model is same as this.rootModel() unless the last element is
    //   this.any or this.every. In which case we have to find
    //   the filter model by digging into the nested arrays.
    let modelName: string = '';
    let anyMap: {[key: string]: boolean} = {};
    let clauses: QueryClause[] = [];

    if (!this.isCollectionModelFilter(rawFilterArray)) {
      modelName = this.rootModel();
      clauses = [
        {
          modelName,
          any: true,
          subclauses: [this.toSubclause(rawFilterArray)]
        }
      ];
    } else {
      // The model name that goes into the QueryFilter has to be
      //   the most nested filter's model (first argument?). This is
      //   typically the value before the innermost array value...
      // Have to also populate the anyMap now, too.
      let lastModelFilter: any[] = [];
      let queue: any[] = [rawFilterArray];
      let activeClauses: {[key: string]: QueryClause} = {};
      let lastActiveClause: QueryClause = {
        modelName: '',
        any: true
      };

      while (queue.length > 0) {
        let filter = queue.shift();

        let firstArg = filter[0];
        if (filter === this.and || filter === this.or) {
          continue;
        } else if (
          Array.isArray(filter) &&
          (firstArg === this.and || firstArg === this.or)
        ) {
          filter.slice(1).forEach((nestedFilter) => {
            queue.push(nestedFilter);
          });
        } else if (this.isCollectionModelFilter(filter)) {
          lastModelFilter = [...filter];
          let currentModelName = lastModelFilter[0]; // is this always true?
          let isAny = lastModelFilter[lastModelFilter.length - 1] === this.any;

          if ('' === modelName) {
            modelName = currentModelName;
            anyMap[currentModelName] = isAny;
          }

          if (!activeClauses[currentModelName]) {
            activeClauses[currentModelName] = {
              modelName: currentModelName,
              any: isAny,
              subclauses: []
            };
          }
          lastActiveClause = activeClauses[currentModelName];
          filter
            .filter((value: any) => Array.isArray(value))
            .forEach((nestedFilter: any[]) => {
              queue.push(nestedFilter);
            });
        } else if (
          QueryFilterOperator.allOperators().some((op) => filter.includes(op))
        ) {
          lastActiveClause.subclauses?.push(this.toSubclause(filter));
        }
      }

      clauses = clauses.concat(Object.values(activeClauses));
    }

    return {
      modelName,
      anyMap,
      clauses
    };
  }

  toSubclause(rawSubclause: any[]): QuerySubclause {
    // For subclauses, I think attributeType will populate itself
    //   in the useQuerySubclauses hook, so we leave it as "" here.
    // For subclauses, some operators have only 2 arguments (::true, ::false, ::untrue). Use this.operatorHelper() to detect these.
    //      - Some of these are even flipped! ::has, ::lacks come before attributeName
    // ::in (and others) requires us to comma-join the array to a single string.
    //
    // rawSubclause must be an array, i.e. ["name", "::=", "foo"]
    //   - It can be nested.
    //   - It can include the path to a different model, instead of just being an attribute on the rootModel()
    //   - It can also start with ::and or ::or ...
    let subclause = {...EmptyQuerySubclause};

    if (this.hasNestedSubclauses(rawSubclause)) {
    } else {
      let testInvertedOperator = rawSubclause[0];
      let testOperator = rawSubclause[1];
      if (
        QueryFilterOperator.terminalInvertOperators.includes(
          testInvertedOperator
        )
      ) {
        subclause.operator = testInvertedOperator;
        subclause.attributeName = rawSubclause[1];
      } else if (QueryFilterOperator.terminalOperators.includes(testOperator)) {
        subclause.operator = testOperator;
        subclause.attributeName = rawSubclause[0];
      } else {
        // Must include an operand!
        subclause.attributeName = rawSubclause[0];
        subclause.operator = rawSubclause[1];

        if (
          QueryFilterOperator.commaSeparatedOperators.includes(
            subclause.operator
          )
        ) {
          subclause.operand = rawSubclause[2].join(',');
        } else {
          subclause.operand = rawSubclause[2];
        }
      }
    }

    return subclause;
  }

  hasNestedSubclauses(rawSubclause: any[]): boolean {
    // Any element except the last one (i.e. from a ::in or ::slice operand)
    //   that is an Array, indicates a nested subclause.
    return rawSubclause
      .slice(0, rawSubclause.length - 1)
      .some((subclause: any) => Array.isArray(subclause));
  }

  orRecordFilterIndices(): number[] {
    // There is no way to get the original order of the filters here, so
    //   we will just have to keep the Or filters at the end of the
    //   filters list. So here we'll just check for the presence of
    //   ::or and return those indices here.
    return this.rawFilters()
      .filter((filter) => this.and !== filter)
      .map((filter, index) => {
        if (Array.isArray(filter) && this.or === filter[0]) {
          return [...Array(filter.length - 1).keys()].map((val) => val + index);
        }
        return null;
      })
      .filter((_) => _)
      .flat() as number[];
  }

  columns(): QueryColumn[] {
    return [];
  }

  rawFilters() {
    // Computed from the query between indices 0 and the location of "::all"
    return this.query.slice(1, this.allIndex()).flat(1);
  }

  rawColumns() {
    // Computed from the query between the location of "::all" and end
    return this.query.slice(this.allIndex() + 1, this.query.length).flat(1);
  }

  allIndex() {
    return this.query.indexOf(this.all);
  }
}
