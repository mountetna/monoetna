import React from 'react';
import {makeStyles} from '@material-ui/core/styles';

import Grid from '@material-ui/core/Grid';
import Typography from '@material-ui/core/Typography';
import Table from '@material-ui/core/Table';
import TableCell from '@material-ui/core/TableCell';
import TableContainer from '@material-ui/core/TableContainer';
import TableHead from '@material-ui/core/TableHead';
import TableRow from '@material-ui/core/TableRow';
import TableBody from '@material-ui/core/TableBody';
import TablePagination from '@material-ui/core/TablePagination';
import AntSwitch from './ant_switch';

import {QueryTableColumn} from '../../contexts/query/query_types';
import QueryTableAttributeViewer from './query_table_attribute_viewer';

const useStyles = makeStyles({
  table: {
    minWidth: 650
  },
  table_container: {
    flex: '1 0',
    overflowY: 'scroll'
  },
  table_controls: {
    padding: '0px 15px'
  },
  hidden_columns: {
    flex: '1 1'
  }
});

const QueryTable = ({
  columns,
  rows,
  numRecords,
  page,
  pageSize,
  maxColumns,
  expandMatrices,
  handlePageChange,
  handlePageSizeChange
}: {
  columns: QueryTableColumn[];
  rows: any;
  numRecords: number;
  page: number;
  pageSize: number;
  maxColumns: number;
  expandMatrices: boolean;
  handlePageChange: (
    e: React.MouseEvent<HTMLButtonElement, MouseEvent> | null,
    newPage: number
  ) => void;
  handlePageSizeChange: (
    e: React.ChangeEvent<HTMLTextAreaElement | HTMLInputElement>
  ) => void;
}) => {
  const classes = useStyles();

  return (
    <React.Fragment>
      <Grid className={classes.table_controls} container>
        <Grid
          item container style={{ width: 'auto', flex: '1 1 auto' }}
          justifyContent='flex-end'
          alignItems='center'
        >
        <Grid className={ classes.hidden_columns } item>
          {columns.length > maxColumns ? (
            <Typography variant='caption' align='right' color='error'>
              {(columns.length - maxColumns).toLocaleString()} columns
              not shown. Add slices to matrix columns or download TSV
              to see entire data frame.
            </Typography>
          ) : null}
        </Grid>
          <TablePagination
            rowsPerPageOptions={[10, 25, 50, 200]}
            component='div'
            count={numRecords}
            rowsPerPage={pageSize}
            page={page}
            onPageChange={handlePageChange}
            onRowsPerPageChange={handlePageSizeChange}
          />
        </Grid>
      </Grid>
      <TableContainer className={classes.table_container}>
        <Table stickyHeader className={classes.table} size='small' aria-label='result table'>
          <TableHead>
            <TableRow>
              {columns
                ?.slice(0, maxColumns)
                .map(({label}: {label: string}, index: number) => (
                  <TableCell key={index} style={{ ...(index == 0) ? { position: 'sticky', left: 0, top: 0, zIndex: 102 } : {} } as React.CSSProperties}>{label}</TableCell>
                ))}
            </TableRow>
          </TableHead>
          <TableBody>
            {rows?.map((row: any[]) => {
              return (
                <TableRow hover tabIndex={-1} key={row[0]}>
                  {row.slice(0, maxColumns).map((datum: any, index: number) => (
                    <TableCell key={index} scope='row' style={{ ...(index == 0) ? { position: 'sticky', left: 0, zIndex: 101, background: 'white'} : {} } as React.CSSProperties}>
                      <QueryTableAttributeViewer
                        tableColumn={columns[index]}
                        expandMatrices={expandMatrices}
                        datum={datum}
                      />
                    </TableCell>
                  ))}
                </TableRow>
              );
            })}
          </TableBody>
        </Table>
      </TableContainer>
    </React.Fragment>
  );
};

export default QueryTable;
