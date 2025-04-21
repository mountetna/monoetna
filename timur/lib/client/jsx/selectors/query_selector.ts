import * as _ from 'lodash';

import {Attribute, Model} from 'etna-js/models/magma-model';
import {
  EmptyQueryClause,
  QueryColumn,
  QuerySlice,
  QueryTableColumn,
  QueryPayload,
  Workflow,
  CreateFigurePayload
} from '../contexts/query/query_types';
import {QueryGraph} from '../utils/query_graph';

export const selectAllowedModelAttributes = (
  attributes: Attribute[],
  includeChildrenModels: boolean = false
): Attribute[] => {
  // Keep "identifier" because it's useful for ::has and ::lacks
  // Don't let folks query "up" the tree, only down it.
  let unallowedAttributeTypes = ['parent', 'link'];

  if (!includeChildrenModels) {
    unallowedAttributeTypes.push('child');
    unallowedAttributeTypes.push('collection');
    unallowedAttributeTypes.push('table');
  }

  return attributes.filter(
    (attr: Attribute) => !unallowedAttributeTypes.includes(attr.attribute_type)
  );
};

export const selectMatrixAttributes = (
  attributes: Attribute[],
  selectedAttributes: QueryColumn[]
): Attribute[] => {
  const selectedAttributeNames = selectedAttributes.map(
    (attr) => attr.attribute_name
  );

  return attributes.filter(
    (attr) =>
      'matrix' === attr.attribute_type &&
      selectedAttributeNames.includes(attr.attribute_name)
  );
};

export const getPath = (
  array: any[],
  heading: string,
  currentPath: number[]
): number[] => {
  if (!array) return [];
  if (!Array.isArray(array)) array = [array];

  let index = array.indexOf(heading);
  if (index > -1) return currentPath.concat([index]);

  let innerPath: number[] = [];
  array.forEach((ele, index: number) => {
    if (Array.isArray(ele)) {
      let tempPath = getPath(ele, heading, currentPath.concat(index));
      if (tempPath.length > 0) {
        innerPath = tempPath;
      }
    }
  });

  return innerPath.length > 0 ? innerPath : [];
};

export const pathToColumn = (
  array: any[],
  heading: string,
  expandMatrices: boolean
): string => {
  let indexlessHeading = heading.split('@')[0];
  let startingIndexPlusMatrixColId = heading.split('@')[1];

  if (!startingIndexPlusMatrixColId) return '-1';

  let fullPath: number[] = [];

  if (expandMatrices) {
    let startingIndex = parseInt(startingIndexPlusMatrixColId.split('.')[0]);
    let sliceColId = startingIndexPlusMatrixColId.split('.')[1];

    fullPath = getPath(array[startingIndex], indexlessHeading, [startingIndex]);

    if (!sliceColId) return fullPath.length > 0 ? fullPath[0].toString() : '-1';

    // fullPath returns the path to project#model::attribute.
    //   The tuple should be [project#model::attribute, [array, of, slice, operands]]
    let pathToSliceOperands = fullPath.slice(0, -1);
    pathToSliceOperands.push(1);
    let sliceOperands = _.at(array, pathToSliceOperands.join('.'))[0];

    if (sliceOperands == null) return '-1';

    let sliceIndex = sliceOperands.indexOf(sliceColId);

    // Get rid of the extra [1] used to find the slice operands
    return pathToSliceOperands.slice(0, -1).concat([sliceIndex]).join('.');
  } else {
    let startingIndex = parseInt(startingIndexPlusMatrixColId);

    if (!Array.isArray(array[startingIndex])) {
      if (array[startingIndex] === indexlessHeading)
        {return startingIndex.toString();}
      else return '-1';
    }

    fullPath = getPath(array[startingIndex], indexlessHeading, [startingIndex]);

    return fullPath.length > 0 ? fullPath[0].toString() : '-1';
  }
};

export const isMatrixSlice = (slice: QuerySlice) => {
  if (!slice.clause.subclauses) return false;

  return slice.clause.subclauses.some((subclause) => {
    return '::slice' === subclause.operator;
  });
};

export const hasMatrixSlice = (column: QueryColumn) => {
  return column.slices.some((slice) => isMatrixSlice(slice));
};

export const emptyQueryClauseStamp = (modelName: string) => {
  return {
    ...EmptyQueryClause,
    modelName
  };
};

export const queryColumnMatrixHeadings = (column: QueryColumn): string[] => {
  return column.slices
    .filter((slice) => isMatrixSlice(slice))
    .map((slice) => {
      return (
        // Given the above isMatrixSlice filter, slice.clause.subclauses
        //   should always exist, but this makes tsc happy.
        (
          (slice.clause.subclauses || []).find((subclause) => {
            return '::slice' === subclause.operator;
          })?.operand as string
        ).split(',')
      );
    })
    .flat()
    .filter((value) => null != value);
};

export const isIdentifierQuery = (
  columns: QueryColumn[] | QueryTableColumn[]
) => {
  return columns.length === 1;
};

export const userColumns = (columns: QueryColumn[]) => {
  const columnLabels = columns.map(
    ({display_label}: {display_label: string}) => display_label
  );

  return isIdentifierQuery(columns)
    ? [columnLabels[0], columnLabels[0]]
    : columnLabels;
};

export const queryPayload = ({
  query,
  columns,
  expandMatrices
}: {
  query: string | any[];
  columns: QueryColumn[];
  expandMatrices: boolean;
}): QueryPayload => {
  return {
    query,
    user_columns: userColumns(columns),
    expand_matrices: expandMatrices,
    format: 'tsv'
  };
};

export const createFigurePayload = ({
  query,
  title,
  workflow
}: {
  query: QueryPayload;
  title: string;
  workflow: Workflow;
}) => {
  let payload: CreateFigurePayload = {
    title,
    workflow_name: workflow.name,
    inputs: {}
  };

  Object.entries(workflow.inputQueryMap || {}).forEach(
    ([cwlInput, inputSource]: [string, string]) => {
      if (!query.hasOwnProperty(inputSource)) return;

      payload.inputs[cwlInput] = query[inputSource as keyof QueryPayload];

      if (
        Array.isArray(payload.inputs[cwlInput]) ||
        (typeof payload.inputs[cwlInput] === 'object' &&
          payload.inputs[cwlInput] !== null)
      ) {
        payload.inputs[cwlInput] = JSON.stringify(payload.inputs[cwlInput]);
      }
    }
  );

  return payload;
};
