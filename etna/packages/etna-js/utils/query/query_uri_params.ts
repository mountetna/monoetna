import {
  QueryFilter,
  QueryClause,
  QueryColumn,
  QuerySlice,
  EmptyQuerySubclause,
  EmptyQueryClause,
  QuerySubclause
} from '../../contexts/query/query_types';

export const isOldClauseFormat = (clause: QueryClause) => {
  return (
    clause.attributeName && clause.attributeType && clause.operator
    // clause.operand could be empty string as a valid entry
  );
};

export const subclauseFromOldClauseFormat = (
  clause: QueryClause
): QuerySubclause => {
  return {
    attributeName: clause.attributeName || '',
    attributeType: clause.attributeType || '',
    operator: clause.operator || '',
    operand: clause.operand || ''
  };
};

export const migrateSubclauses = (
  recordFilters: QueryFilter[]
): QueryFilter[] => {
  return recordFilters.map((filter: QueryFilter) => {
    return {
      ...filter,
      clauses: filter.clauses.map((clause: QueryClause) => {
        let subclauses: QuerySubclause[] = [];
        if (isOldClauseFormat(clause)) {
          subclauses.push(subclauseFromOldClauseFormat(clause));
          subclauses = subclauses.concat(clause.subclauses || []);
        } else {
          subclauses = [
            ...(clause.subclauses || [
              {
                ...EmptyQuerySubclause
              }
            ])
          ];
        }

        return {
          modelName: clause.modelName,
          any: clause.any,
          subclauses
        };
      })
    };
  });
};

export const migrateSlices = (columns: QueryColumn[]): QueryColumn[] => {
  return columns.map((column: QueryColumn) => {
    return {
      ...column,
      slices: column.slices.map((slice: QuerySlice) => {
        let clause: QueryClause = {
          modelName: slice.clause.modelName,
          any: slice.clause.any
        };
        if (isOldClauseFormat(slice.clause)) {
          clause.subclauses = [subclauseFromOldClauseFormat(slice.clause)];
        } else {
          clause.subclauses = [
            ...(slice.clause.subclauses || [
              {
                ...EmptyQuerySubclause
              }
            ])
          ];
        }

        return {
          ...slice,
          clause
        };
      })
    };
  });
};
