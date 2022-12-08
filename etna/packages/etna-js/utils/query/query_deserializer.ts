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

import {range} from 'lodash';

export class QueryDeserializer {
  query: any[];
  userColumns: string[];

  constructor(query: any[], userColumns: string[]) {
    this.query = query;
    this.userColumns = userColumns;
  }

  rootModel() {
    return this.query[0];
  }

  recordFilters() {}

  orRecordFilterIndices() {
    // There is no way to get the original order of the filters here, so
    //   we will just have to keep the Or filters at the end of the
    //   filters list. So here we'll just check for the presence of
    //   ::or and return those indices here.
    this.rawFilters()
      .filter((filter) => filter !== '::and')
      .map((filter, index) => {
        if (Array.isArray(filter) && filter[0] === '::or') {
          return [...Array(filter.length).keys()].map((val) => val + index);
        }
        return null;
      })
      .filter((_) => _)
      .flat();
  }

  columns() {}

  rawFilters() {
    // Computed from the query between indices 0 and the location of "::all"
    return this.query.slice(1, this.query.indexOf(':all'));
  }

  rawColumns() {
    // Computed from the query between the location of "::all" and end
    return this.query.slice(this.query.indexOf('::all'), this.query.length);
  }
}
