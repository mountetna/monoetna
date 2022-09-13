import {
  QueryFilter,
  QueryClause,
  EmptyQuerySubclause
} from '../contexts/query/query_types';

export const migrateSubclauses = (
  recordFilters: QueryFilter[]
): QueryFilter[] => {
  return recordFilters.map((filter: QueryFilter) => {
    return {
      ...filter,
      clauses: filter.clauses.map((clause: QueryClause) => {
        let subclauses = [];
        if (
          clause.attributeName &&
          clause.attributeType &&
          clause.operator &&
          clause.operand
        ) {
          subclauses.push({
            attributeName: clause.attributeName,
            attributeType: clause.attributeType,
            operator: clause.operator,
            operand: clause.operand
          });
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
