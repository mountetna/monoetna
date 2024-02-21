import React, {useState, useCallback, useMemo, useEffect} from 'react';
import Typography from '@material-ui/core/Typography';
import LinearProgress from '@material-ui/core/LinearProgress';
import TextField from '@material-ui/core/TextField';
import Table from '@material-ui/core/Table';
import TableBody from '@material-ui/core/TableBody';
import TableCell from '@material-ui/core/TableCell';
import TableSortLabel from '@material-ui/core/TableSortLabel';
import TableHead from '@material-ui/core/TableHead';
import TableRow from '@material-ui/core/TableRow';
import TableContainer from '@material-ui/core/TableContainer';
import Checkbox from '@material-ui/core/Checkbox';
import Grid from '@material-ui/core/Grid';
import {makeStyles} from '@material-ui/core/styles';
import {isEqual} from 'lodash';
import {sortAttributeList} from '../../utils/attributes';
import VisibilityOffIcon from '@material-ui/icons/VisibilityOff';
import Tooltip from '@material-ui/core/Tooltip';
import InputAdornment from '@material-ui/core/InputAdornment';
import SearchIcon from '@material-ui/icons/Search';
import Paper from '@material-ui/core/Paper';
import {Template,Attribute} from '../../api/magma_api';

const attributeStyles = makeStyles((theme) => ({
  attribute_row: {
    height: '36px'
  },
  filter: {
    paddingBottom: '10px'
  },
  cell: {
  },
  value: {
    color: 'darkgoldenrod',
    cursor: 'pointer'
  },
  missing: {},
  indicator: {
    width: '30px',
    color: 'gray',
    cursor: 'default'
  },
  add: {
    width: '30px',
    color: 'green'
  },
  type: {
    color: 'gray',
    width: '120px',
    paddingRight: '10px'
  },
  progress: {
    width: '30px',
    margin: '5px'
  },
  description: {
    width: '300px'
  },
  descbox: {
    width: '300px',
    height: '25px',
    overflow: 'hidden',
    whiteSpace: 'nowrap',
    textOverflow: 'ellipsis'
  },
  counts: {
    minWidth: '120px'
  },
  ident: {},
  changed: {
    backgroundColor: 'rgba(255,255,0,0.1)'
  },
  present: {
    backgroundColor: 'rgba(0,255,0,0.1)'
  },
  absent: {
    backgroundColor: 'rgba(255,0,0,0.1)'
  },
  hiddenIcon: {
    marginRight: '1rem'
  },
  typeWrapper: {
    display: 'flex',
    justifyContent: 'right'
  },
  table: {
    boxSizing: 'border-box',
    height: 'auto'
  },
  hiddenTypeWrapper: {
    display: 'flex',
    justifyContent: 'right',
    paddingTop: '0.35rem'
  }
}));

const diffTypes = {
  ident: {ind: '', title: ''},
  present: {ind: '+', title: 'Present in this model'},
  absent: {ind: '-', title: 'Absent in this model'},
  changed: {ind: 'c', title: 'Changed in this model'}
};

const ATT_KEYS:{ [key: string]: keyof Attribute } = {
  attribute: 'attribute_name',
  type: 'attribute_type',
  group: 'attribute_group',
  description: 'description'
};

const ModelAttribute = ({
  attribute_name,
  template,
  diffTemplate,
  setAttribute,
  count,
  modelCount,
  selectAttribute,
  selected,
  isHidden
}:{
  attribute_name: string;
  template: Template|null;
  diffTemplate?: Template;
  setAttribute?: Function;
  count?: number;
  modelCount?: number;
  selectAttribute?: Function;
  selected: boolean;
  isHidden: boolean;
}) => {
  if (!template) return null;

  const classes = attributeStyles();

  const attribute = template.attributes[attribute_name];

  const diffAttribute = diffTemplate?.attributes[attribute_name];

  const displayAttribute = diffTemplate && !attribute ? diffAttribute : attribute;

  if (!displayAttribute) return null;

  const diffType = !diffTemplate
    ? 'ident'
    : attribute && diffAttribute
      ? isEqual(attribute, diffAttribute) ? 'ident' : 'changed'
      : attribute ? 'present' : 'absent';

  const {attribute_type, attribute_group, description} = displayAttribute;

  return (
    <TableRow className={`${classes.attribute_row} ${classes[diffType]}`}>
      {selectAttribute &&
        <TableCell padding="checkbox">
          <Checkbox
            checked={selected}
            onChange={() => selectAttribute(attribute_name)}
          />
        </TableCell>
      }
      {diffTemplate && (
        <TableCell
          className={classes.indicator}
          align='left'
          title={diffTypes[diffType].title}
        >
          {diffTypes[diffType].ind}
        </TableCell>
      )}
      <TableCell className={classes.type} align='right'>
        <div
          className={isHidden ? classes.hiddenTypeWrapper : classes.typeWrapper}
        >
          {isHidden ? (
            <div className={classes.hiddenIcon}>
              <VisibilityOffIcon fontSize='small'/>
            </div>
          ) : null}
          {attribute_type}
        </div>
      </TableCell>
      <TableCell
        className={attribute ? classes.value : classes.missing}
        align='left'
        onClick={(attribute && setAttribute) ? () => setAttribute(attribute_name) : undefined}
      >
        {attribute_name}{' '}
      </TableCell>
      <TableCell align='left'>{attribute_group}</TableCell>
      <TableCell className={classes.description} align='left'>
        <Tooltip title={description || ''} aria-label='Description'>
          <div className={classes.descbox}>{description}</div>
        </Tooltip>
      </TableCell>
      {count != undefined && (
        <TableCell className={classes.counts} align='left'>
          <Grid container alignItems='center'>
            {count}
            <LinearProgress
              className={classes.progress}
              variant='determinate'
              value={(100 * count) / (modelCount || 1)}
            />
            <Typography variant='body2' color='textSecondary'>
              {Math.round((100 * count) / (modelCount || 1))}%
            </Typography>
          </Grid>
        </TableCell>
      )}
    </TableRow>
  );
};

const ModelAttributesTable = ({
  template, diffTemplate, attributeCounts, showHiddenAttributes,
  setAttribute, modelCount, selected={}, setSelected, className
}:{
  template: Template|null;
  diffTemplate?: Template;
  attributeCounts?: { [key: string]: number };
  showHiddenAttributes?: boolean;
  setAttribute?: Function;
  count?: number;
  modelCount?: number;
  selected: { [key: string]: boolean };
  setSelected?: Function;
  className?: string;
}) => {
  const classes = attributeStyles();
  const [order, setOrder] = useState<'asc'|'desc'>('asc');
  const [orderBy, setOrderBy] = useState('type');
  const [filterString, setFilterString] = useState('');

  const sortByOrder = (attributes:Attribute[]) => {
    let srt;
    if (orderBy === 'type') srt = sortAttributeList(attributes);
    else {
      srt = attributes.sort((a, b) =>
        ((a[ATT_KEYS[orderBy]] || '') as string).localeCompare((b[ATT_KEYS[orderBy]] || '') as string)
      );
    }
    return order === 'desc' ? srt.reverse() : srt;
  };

  const attributes = Object.values({
    ...(template?.attributes || {}),
    ...(diffTemplate && diffTemplate.attributes)
  }).filter(
    attribute => showHiddenAttributes || !attribute.hidden
  );

  const filterMatch = new RegExp(
    `^(?:(${Object.keys(ATT_KEYS).join('|')}):)?(.*)$`
  );

  const matchesFilter = useCallback(
    (attribute) => {
      if (filterString == '') return true;

      let tokens = filterString
        .split(/\s/)
        .filter(_ => _)
        .map((token) => token.match(filterMatch)?.slice(1) as string[]);

      return tokens.every(([column, token]) => {
        const tokenMatch = new RegExp(token, 'i');
        const values = column
          ? [attribute[ATT_KEYS[column]]]
          : Object.values(ATT_KEYS).map((k) => attribute[k]);

        return values.some((s) => s?.match(tokenMatch));
      });
    },
    [filterString]
  );

  const selectAttribute = (attribute_name:string) => {
    if (selected[attribute_name]) {
      const { [attribute_name]: _, ...newSelected } = selected;
      (setSelected as Function)(newSelected);
    } else
      (setSelected as Function)({ [attribute_name]: true, ...selected } );
  }

  const numSelected = Object.keys(selected).length;

  const visibleSelected = attributes.filter(a => selected[a.attribute_name] && matchesFilter(a));

  const visibleRowCount = attributes.filter(matchesFilter).length;

  const selectAll = () => {
    const newSelected = {...selected};
    if (visibleSelected.length == visibleRowCount) {
      // deselect all visible
      visibleSelected.forEach(a => delete newSelected[a.attribute_name]);
    } else {
      // select all visible
      attributes.filter(
        a => !selected[a.attribute_name] && matchesFilter(a)
      ).forEach(
        a => newSelected[a.attribute_name] = true
      );
    }
    (setSelected as Function)(newSelected);
  }

  return <Grid className={className}>
    <TextField
      fullWidth
      placeholder='Filter attributes, e.g. "rna type:file"'
      variant='outlined'
      size='small'
      className={classes.filter}
      value={filterString}
      InputProps={{
        startAdornment: (
          <InputAdornment position='start'>
            <SearchIcon />
          </InputAdornment>
        )
      }}
      onChange={(e) => setFilterString(e.target.value)}
    />
    <TableContainer component={Paper} className={classes.table}>
      <Table stickyHeader size='small'>
        <TableHead>
          <TableRow>
            {setSelected &&
              <TableCell padding="checkbox">
                <Checkbox
                  indeterminate={visibleSelected.length > 0 && visibleSelected.length < visibleRowCount}
                  checked={visibleRowCount > 0 && visibleSelected.length === visibleRowCount}
                  onChange={selectAll}
                  inputProps={{ 'aria-label': 'select all desserts' }}
                />
              </TableCell>
            }
            {diffTemplate && <TableCell className={classes.indicator} />}
            {['type', 'attribute', 'group', 'description', 'counts'].map(
              (column_name) =>
                (column_name !== 'counts' || attributeCounts) && (
                  <TableCell
                    key={column_name}
                    className={column_name in classes ? classes[column_name as keyof typeof classes] : undefined}
                    align={column_name == 'type' ? 'right' : 'left'}
                  >
                    <TableSortLabel
                      active={orderBy === column_name}
                      direction={orderBy === column_name ? order : 'asc'}
                      onClick={() => {
                        setOrder(
                          orderBy === column_name && order === 'asc' ? 'desc' : 'asc'
                        );
                        setOrderBy(column_name);
                      }}
                    >
                      {column_name}
                    </TableSortLabel>
                  </TableCell>
                )
            )}
          </TableRow>
        </TableHead>

        <TableBody>
          {sortByOrder(attributes).filter(matchesFilter)
            .map((attribute:Attribute) => 
              <ModelAttribute
                key={attribute.attribute_name}
                attribute_name={attribute.attribute_name}
                setAttribute={setAttribute}
                count={
                  attributeCounts && attributeCounts[attribute.attribute_name]
                }
                modelCount={modelCount}
                template={template}
                diffTemplate={diffTemplate}
                isHidden={!!attribute.hidden}
                selected={selected[attribute.attribute_name] || false}
                selectAttribute={setSelected ? selectAttribute : undefined}
              />
            )}
        </TableBody>
      </Table>
    </TableContainer>
  </Grid>
}

export default ModelAttributesTable
