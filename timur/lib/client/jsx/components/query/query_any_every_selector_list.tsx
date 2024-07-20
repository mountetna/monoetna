import React, {useState, useCallback} from 'react';

import Selector from './query_selector';
import MenuItem from '@material-ui/core/MenuItem';
import Typography from '@material-ui/core/Typography';

import {makeStyles} from '@material-ui/core/styles';

import {QueryFilter} from '../../contexts/query/query_types';
import {useEffect} from 'react';

const useStyles = makeStyles((theme) => ({
  textInput: {
    margin: theme.spacing(1),
    minWidth: 120,
    paddingLeft: '1rem'
  }
}));

const QueryAnyEverySelectorList = ({
  filter,
  index,
  patchRecordFilter
}: {
  filter: QueryFilter;
  index: number;
  patchRecordFilter: (index: number, updatedFilter: QueryFilter) => void;
}) => {
  const [anyMap, setAnyMap] = useState({} as {[key: string]: boolean});
  const classes = useStyles();

  const handlePatchFilter = useCallback(
    (modelName: string) => {
      let copy = {
        ...filter,
        anyMap: {
          ...filter.anyMap,
          [modelName]: !filter.anyMap[modelName]
        }
      };
      patchRecordFilter(index, copy);
    },
    [patchRecordFilter, index, filter]
  );

  useEffect(() => {
    setAnyMap(filter.anyMap);
  }, [filter.anyMap]);

  const anyList = Object.entries(anyMap || {});

  return (
    <React.Fragment>
      <Typography>{anyList.length ? 'For' : 'For the'}&nbsp;</Typography>
      {anyList.map(([modelName, value], index: number) => {
        return ( <>
          <Selector
            canEdit={true}
            key={index}
            name={value.toString()}
            onSelect={() => handlePatchFilter(modelName)}
            choiceSet={ [
              [ `any ${ (index < anyList.length - 1) ? modelName : ''}`, true ],
              [ `every ${ (index < anyList.length - 1) ? modelName : ''}`, false ]
            ]}
          />
          {
            (index < anyList.length - 1) && <Typography>&nbsp;and&nbsp;</Typography>
          }
        </>
        );
      })}
    </React.Fragment>
  );
};

export default QueryAnyEverySelectorList;
