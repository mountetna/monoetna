import React, {useCallback} from 'react';

import downloadjs from 'downloadjs';

import {useModal} from 'etna-js/components/ModalDialogContainer';
import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {showMessages} from 'etna-js/actions/message_actions';
import {requestAnswer} from 'etna-js/actions/magma_actions';
import {
  QueryResponse,
  EmptyQueryResponse,
  QueryColumn
} from '../../contexts/query/query_types';
import {Cancellable} from 'etna-js/utils/cancellable';
import {queryPayload} from '../../selectors/query_selector';

const useResultsActions = ({
  countQuery,
  query,
  page,
  pageSize,
  columns,
  expandMatrices,
  showDisconnected,
  setDataAndNumRecords
}: {
  countQuery: string | any[];
  query: string | any[];
  page: number;
  pageSize: number;
  columns: QueryColumn[];
  expandMatrices: boolean;
  showDisconnected: boolean;
  setDataAndNumRecords: (data: QueryResponse, count: number) => void;
}) => {
  const invoke = useActionInvoker();
  const {dismissModal} = useModal();

  const runQuery = useCallback((queryPage?: number) => {
    if ('' === countQuery || '' === query) return;

    let numRecords = 0;
    
    if (queryPage === undefined) queryPage = page;

    setDataAndNumRecords(EmptyQueryResponse, 0);
    invoke(requestAnswer({query: countQuery, show_disconnected: showDisconnected}))
      .then((countData: {answer: number}) => {
        numRecords = countData.answer;
        return invoke(
          requestAnswer({query, page_size: pageSize, page: queryPage + 1, show_disconnected: showDisconnected})
        );
      })
      .then((answerData: QueryResponse) => {
        setDataAndNumRecords(answerData, numRecords);
      })
      .catch((e: any) => {
        Promise.resolve(e).then((error: {[key: string]: string[]}) => {
          console.error(error);
          if (error.errors?.[0] == 'No results found') {
            setDataAndNumRecords(EmptyQueryResponse, 0);
          }
          invoke(showMessages(error.errors || [error.error] || error));
        });
      });
  }, [query, countQuery, pageSize, page, invoke, setDataAndNumRecords]);

  const downloadData = useCallback(
    ({transpose}: {transpose: boolean}) => {
      if ('' === query) return;

      const cancellable = new Cancellable();

      cancellable
        .race(
          invoke(
            requestAnswer({
              ...queryPayload({query, columns, expandMatrices}), show_disconnected: showDisconnected,
              transpose
            })
          )
        )
        .then(({result, cancelled}: any) => {
          if (result && !cancelled) {
            downloadjs(
              result,
              `${
                CONFIG.project_name
              }-query-results-${new Date().toISOString()}.tsv`,
              'text/tsv'
            );
            dismissModal();
          }
        })
        .catch((error: any) => {
          Promise.resolve(error).then((e) => {
            console.error(e);
            invoke(showMessages(e.errors || [e.toString()]));
          });
        });

      return () => cancellable.cancel();
    },
    [query, columns, invoke, expandMatrices, dismissModal]
  );

  return {
    runQuery,
    downloadData
  };
};

export default useResultsActions;
