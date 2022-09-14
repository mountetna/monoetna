import {
  QueryFilter,
  QueryClause,
  EmptyQuerySubclause
} from '../contexts/query/query_types';

export const isOldClauseFormat = (clause: QueryClause) => {
  return (
    clause.attributeName && clause.attributeType && clause.operator
    // clause.operand could be empty string as a valid entry
  );
};

export const cloneOldClauseFormat = (clause: QueryClause) => {
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
        let subclauses = [];
        if (isOldClauseFormat(clause)) {
          subclauses.push(cloneOldClauseFormat(clause));
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
