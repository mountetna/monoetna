import {
  QueryFilter,
  QueryClause,
  QueryColumn,
  QuerySlice,
  EmptyQuerySubclause,
  EmptyQueryClause,
  QuerySubclause
} from '../contexts/query/query_types';

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

const btou = (str: string) => btoa(str).replace(/\+/g, '-').replace(/\//g, '_').replace(/=/g, '');
const utob = (str: string) => {
  str = str.replace(/-/g, '+').replace(/_/g, '/');
  while (str.length % 4) str += '=';
  return atob(str);
}

export const decodeCompressedParams = async (params: string) => {
  const stream = new Blob([
    Uint8Array.from(
      utob(params),
      c => c.charCodeAt(0)
    )
  ]).stream();
  let decompressedStream = stream.pipeThrough(
    new DecompressionStream("gzip")
  );
  const blob = await new Response(decompressedStream as BodyInit).blob();
  const text = await blob.text();
  return JSON.parse(text);
};

export const unpackParams = async (search: string) => {
  let searchParams = new URLSearchParams(search);

  if (searchParams.has('q')) {
    return await decodeCompressedParams(searchParams.get('q') as string);
  }

  return {
    rootModel: searchParams.get('rootModel') || '',
    recordFilters: migrateSubclauses(
      JSON.parse(searchParams.get('recordFilters') || '[]')
    ),
    orRecordFilterIndices: JSON.parse(
      searchParams.get('orRecordFilterIndices') || '[]'
    ),
    columns: migrateSlices(JSON.parse(searchParams.get('columns') || '[]'))
  }
};

export const packParams = async (params: any) => {
  const stream = new Blob(
    [JSON.stringify(params)],
    { type: 'applicaton/json' }
  ).stream();

  const gzipStream = stream.pipeThrough(new CompressionStream("gzip"));

  const compressedResponse = await new Response(gzipStream as BodyInit);
  const blob = await compressedResponse.blob();
  const buffer = await blob.arrayBuffer();
  const compressedBase64 = btou(
    String.fromCharCode(
      ...new Uint8Array(buffer)
    )
  );

  const searchParams = new URLSearchParams();

  searchParams.set('q', compressedBase64)

  return searchParams.toString();
}
