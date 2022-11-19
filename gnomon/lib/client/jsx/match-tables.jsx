import React, {useState, useEffect } from 'react';
import Grid from '@material-ui/core/Grid';
import Link from '@material-ui/core/Link';
import Typography from '@material-ui/core/Typography';
import ProjectHeader from 'etna-js/components/project-header';
import {makeStyles} from '@material-ui/core/styles';
import Letter from './letter';
import Bracket from './bracket';
import Corner from './corner';
import { json_get, json_post } from 'etna-js/utils/fetch';
import { dateFormat, authorFormat } from 'etna-js/utils/format';
import { magmaPath } from 'etna-js/api/magma_api';
import {getDocuments} from 'etna-js/api/magma_api';
import TokenEditor from './token-editor';
import CheckBoxOutlinedIcon from '@material-ui/icons/CheckBoxOutlined';
import LaunchIcon from '@material-ui/icons/Launch';
import Table from '@material-ui/core/Table';
import TableBody from '@material-ui/core/TableBody';
import TableCell from '@material-ui/core/TableCell';
import TableSortLabel from '@material-ui/core/TableSortLabel';
import TableContainer from '@material-ui/core/TableContainer';
import TableHead from '@material-ui/core/TableHead';
import TableRow from '@material-ui/core/TableRow';
import Paper from '@material-ui/core/Paper';
import Tooltip from '@material-ui/core/Tooltip';
import Toolbar from '@material-ui/core/Toolbar';
import Button from '@material-ui/core/Button';

const useStyles = makeStyles(theme => ({
  highlight: {
    background: '#666',
    '& .MuiTableCell-body': {
      color: 'white'
    }
  },
  highlight2: {
    background: '#ccc',
    '& .MuiTableCell-body': {
      color: 'black'
    }
  },
  table_header: {
    borderBottom: '1px solid #ccc',
    minHeight: 'auto',
    padding: '0px 15px',
    background: '#eee'
  }
}))

export const TableWithTitle = ({title, className, children}) => {
  const classes = useStyles();
  return <Grid className={className} item>
    <Toolbar disableGutters={true} className={classes.table_header}>
      <Typography variant='h6' >{title}</Typography>
    </Toolbar>
    { children }
  </Grid>
}

const Rule = ({rule, rule_name, project_name, markNotCreated, highlight, highlight2}) => {
  const classes = useStyles();
  return <TableRow className={(highlight2 && !rule.name_created_at) ? classes.highlight2 : highlight ? classes.highlight : ''}>
    <TableCell component='th' scope='row'>
      {rule_name}
    </TableCell>
    <TableCell>{rule.name + (!markNotCreated || rule.name_created_at ? '' : '*')}</TableCell>
    <TableCell align='center'>{dateFormat(rule.name_created_at, null)}</TableCell>
    <TableCell align='center'>
      {
        rule.record_created_at
        ?  <Link
              color='secondary'
              href={`${CONFIG.timur_host}/${project_name}/browse/${rule_name}/${rule.name}`}>
              {dateFormat(rule.record_created_at)}
            </Link>
          : ''
      }
    </TableCell>
  </TableRow>
}

const getIdentifier = decomposition => decomposition ? decomposition.tokens.map(t => t[1]).join('') : null;

export const IdTreeTable = ({decomposition, project_name, markNotCreated = false, highlight=false, ...props}) => {
  const primary_rule = decomposition != null ? decomposition.rule_name : null;
  const identifier = getIdentifier(decomposition);
  const classes = useStyles();

  return <MatchTable
    headCells={[ 
      {
        name: 'Rule',
        key: 'rule_name',
        contents: rule => rule.rule_name
      },
      {
        name: 'Identifier',
        key: 'name',
        contents: rule => rule.name + (!markNotCreated || rule.name_created_at ? '' : '*')
      },
      {
        name: 'Created',
        align: 'center',
        title: 'Identifiers previously assigned for a current or planned entity',
        contents: rule => dateFormat(rule.name_created_at, null)
      },
      {
        name: 'Recorded',
        align: 'center',
        title: 'Identifiers with data in the Data Library',
        contents: rule => rule.record_created_at ? <Link
          color='secondary'
          href={`${CONFIG.timur_host}/${project_name}/browse/${rule.rule_name}/${rule.name}`}>
          {dateFormat(rule.record_created_at)}
        </Link> : ''
      }
    ]}
    items={!decomposition ? [] : Object.keys(decomposition.rules).map( rule_name => ({ ...decomposition.rules[rule_name], rule_name}))}
    itemClassName={ rule => (highlight && !rule.name_created_at) ? classes.highlight2 : rule.name == identifier ? classes.highlight : ''}
    emptyText='Identifier is incomplete'
    csvName={ `${primary_rule}-matching-rules-${(new Date()).toISOString().slice(0,10)}.csv` }
    {...props}
  />
  //decomposition != null ? Object.keys(decomposition.rules).sort( (a,b) => b==primary_rule ? 1 : a==primary_rule ? -1 : 0 ).map(
}

const Author = ({author}) => {
  const { name, email } = authorFormat(author);
  return <Tooltip title={email}>
    <>{name || email}</>
  </Tooltip>
}

export const MatchingNamesTable = ({names, rule_name, decomposition, ...props}) => {
  const classes = useStyles();
  if (names == null) return null
  const identifier = getIdentifier(decomposition);

  return <MatchTable
    headCells={ [
      {
        name: 'Identifier',
        key: 'identifier',
        contents: name => name.identifier
      },
      {
        name: 'Author',
        align: 'left',
        key: 'author',
        contents: name => <Author author={name.author}/>
      },
      {
        name: 'Created',
        align: 'right',
        key: 'name_created_at',
        contents: name => dateFormat(name.name_created_at, null)
      },
      {
        name: 'Recorded',
        align: 'right',
        key: 'record_created_at',
        contents: name => dateFormat(name.record_created_at, null)
      }
    ] }
    items={names}
    itemClassName={ name => identifier == name.identifier ? classes.highlight : '' }
    emptyText='None'
    csvName={ `${rule_name}-matching-names-${(new Date()).toISOString().slice(0,10)}.csv` }
    {...props}
  />
}

export const MatchTable = ({headCells, csvName, tableClassName, boxClassName, title, items, itemClassName, emptyText }) => {
  const [ orderBy, setOrderBy ] = useState(headCells[0].key);
  const [ order, setOrder ] = useState('asc');
  const classes = useStyles();

  const show = items && items.length > 0;

  const download = () => {
    const data = items.map(
      item => headCells.map( head => item[ head.key ] ).join(',')
    ).join('\n');
    const blob = new Blob([data], { type: 'text/csv' });
    const a = document.createElement('a')
    a.setAttribute('href', window.URL.createObjectURL(blob) )
    a.setAttribute('download', csvName);
    a.click()
  }

  return <Grid className={boxClassName} item>
    <Toolbar disableGutters={true} className={classes.table_header}>
      <Grid item xs={10}><Typography variant='h6' >{title}</Typography></Grid>
      { show && <Grid item align='right' xs={2}><Button onClick={ download }>TSV</Button></Grid> }
    </Toolbar>
    <TableContainer className={tableClassName}>
      <Table stickyHeader size='small'>
        <TableHead>
          {
             show && <TableRow>
              {
                headCells.map(
                  head => <TableCell key={head.key} align={head.align}>
                    <TableSortLabel
                      active={orderBy === head.key}
                      direction={orderBy === head.key ? order : 'asc'}
                      onClick={(key => e => {
                        setOrderBy(key);
                        setOrder( orderBy === key && order === 'asc' ? 'desc' : 'asc' );
                      })(head.key)}>
                    {
                      head.tooltip ?
                        <Tooltip title={head.title} placement='top'>
                          <Typography component='span'>{head.name}</Typography>
                        </Tooltip>
                        : <Typography component='span'>{head.name}</Typography>
                    }
                    </TableSortLabel>
                  </TableCell>
                )
              }
            </TableRow>
          }
        </TableHead>
        {
          !show ?
            <TableBody>
              <TableRow>
                <TableCell colSpan={5}><Typography variant='body1' color='secondary' component='em'>{emptyText}</Typography></TableCell>
              </TableRow>
            </TableBody>
            : <TableBody>
              {items.sort( (i1, i2) => (order == 'desc' ? 1 : -1) * (i1[ orderBy ]||'').localeCompare((i2[orderBy]||''))).map((item, i) => (
              <TableRow key={i} className={itemClassName(item)}>
                {
                  headCells.map(
                    head => <TableCell key={head.key} align={head.align}>{ head.contents( item ) }</TableCell>
                  )
                }
              </TableRow>
            ))}
          </TableBody>
        }
      </Table>
    </TableContainer>
  </Grid>
}
