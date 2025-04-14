import React, {useMemo, useContext} from 'react';

import AttributeViewer from '../attributes/attribute_viewer';
import {QueryGraph} from '../../utils/query_graph';
import {QueryGraphContext} from '../../contexts/query/query_graph_context';
import {QueryTableColumn} from '../../contexts/query/query_types';

const QueryTableAttributeViewer = ({
  tableColumn,
  datum,
  expandMatrices
}: {
  tableColumn: QueryTableColumn;
  datum: any;
  expandMatrices: boolean;
}) => {
  function filename(path: string | null) {
    return path == null
      ? null
      : new URL(`https://${path}`).pathname.split('/').pop();
  }

  const { state: {graph} } = useContext(QueryGraphContext);

  const {attribute, modelName, matrixHeadings, predicate} = tableColumn;

  const isMd5Predicate = useMemo(() => predicate && [ 'size', 'md5' ].includes(predicate), [predicate]);

  const mockRecord = useMemo(() => {
    const {attribute} = tableColumn;

    if (!attribute) return {};

    switch (attribute.attribute_type) {
      case 'file':
      case 'image':
        let name = filename(datum);

        return {
          [attribute.attribute_name]: isMd5Predicate
            ? datum
            : name
            ? {
                url: datum,
                original_filename: name,
                path: name
              }
            : null
        };
      case 'file_collection':
        return {
          [attribute.attribute_name]: isMd5Predicate
            ? datum
            : datum?.map((datum: string) => {
                return {
                  url: datum,
                  original_filename: filename(datum),
                  path: filename(datum)
                };
              })
        };
      default:
        return {
          [attribute.attribute_name]: datum
        };
    }
  }, [tableColumn, datum, isMd5Predicate]);

  if (!tableColumn.attribute) return null;

  return (
    <>
      {expandMatrices && attribute.attribute_type == 'matrix' ? (
        datum?.toString()
      ) : (
        <AttributeViewer
          attribute_name={attribute.attribute_name}
          record={mockRecord}
          model_name={modelName}
          template={graph.template(modelName)}
          mode='model_viewer'
          sliceValues={
            matrixHeadings && matrixHeadings.length > 0 ? matrixHeadings : null
          }
          predicate={predicate}
        />
      )}
    </>
  );
};

export default QueryTableAttributeViewer;
