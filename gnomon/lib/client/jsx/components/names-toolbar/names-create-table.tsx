import React, { useState } from "react";
import { batch } from 'react-redux';
import { createStyles, makeStyles, Theme } from '@material-ui/core/styles';
import Table from '@material-ui/core/Table';
import TableBody from '@material-ui/core/TableBody';
import TableCell from '@material-ui/core/TableCell';
import TableContainer from '@material-ui/core/TableContainer';
import TableHead from '@material-ui/core/TableHead';
import TablePagination from '@material-ui/core/TablePagination';
import TableRow from '@material-ui/core/TableRow';
import TableSortLabel from '@material-ui/core/TableSortLabel';



// adapted from https://v4.mui.com/components/tables/#sorting-amp-selecting
// switch to @mui/x-data-grid when react >= v17

export interface Data {
    name: string
    rule: string
    implicit: boolean
}

function descendingComparator<T>(a: T, b: T, orderBy: keyof T) {
    if (b[orderBy] < a[orderBy]) {
        return -1;
    }
    if (b[orderBy] > a[orderBy]) {
        return 1;
    }
    return 0;
}

type Order = 'asc' | 'desc';

function getComparator<Key extends keyof any>(
    order: Order,
    orderBy: Key,
): (a: { [key in Key]: number | string }, b: { [key in Key]: number | string }) => number {
    return order === 'desc'
        ? (a, b) => descendingComparator(a, b, orderBy)
        : (a, b) => -descendingComparator(a, b, orderBy);
}

function stableSort<T>(array: T[], comparator: (a: T, b: T) => number) {
    const stabilizedThis = array.map((el, index) => [el, index] as [T, number]);
    stabilizedThis.sort((a, b) => {
        const order = comparator(a[0], b[0]);
        if (order !== 0) return order;
        return a[1] - b[1];
    });
    return stabilizedThis.map((el) => el[0]);
}

interface HeadCell {
    id: keyof Data
    label: string
}

const headCells: HeadCell[] = [
    { id: 'name', label: 'Name' },
    { id: 'rule', label: 'Rule' },
];

interface EnhancedTableProps {
    classes: ReturnType<typeof useStyles>;
    onRequestSort: (event: React.MouseEvent<unknown>, property: keyof Data) => void;
    order: Order;
    orderBy: string;
}

function EnhancedTableHead(props: EnhancedTableProps) {
    const { classes, order, orderBy, onRequestSort } = props;
    const createSortHandler = (property: keyof Data) => (event: React.MouseEvent<unknown>) => {
        onRequestSort(event, property);
    };

    return (
        <TableHead>
            <TableRow>
                {headCells.map((headCell) => (
                    <TableCell
                        key={headCell.id}
                        align={'left'}
                        padding={'none'}
                        sortDirection={orderBy === headCell.id ? order : false}
                    >
                        <TableSortLabel
                            active={orderBy === headCell.id}
                            direction={orderBy === headCell.id ? order : 'asc'}
                            onClick={createSortHandler(headCell.id)}
                        >
                            {headCell.label}
                            {orderBy === headCell.id ? (
                                <span className={classes.visuallyHidden}>
                                    {order === 'desc' ? 'sorted descending' : 'sorted ascending'}
                                </span>
                            ) : null}
                        </TableSortLabel>
                    </TableCell>
                ))}
            </TableRow>
        </TableHead>
    );
}

const useStyles = makeStyles((theme: Theme) =>
    createStyles({
        root: {
            width: '100%',
        },
        paper: {
            width: '100%',
            marginBottom: theme.spacing(2),
        },
        table: {
            // minWidth: 750,
        },
        visuallyHidden: {
            border: 0,
            clip: 'rect(0 0 0 0)',
            height: 1,
            margin: -1,
            overflow: 'hidden',
            padding: 0,
            position: 'absolute',
            top: 20,
            width: 1,
        },
        implicitData: {
            fontStyle: "italic",
        },
    }),
);


const NamesTable = ({ rows }: { rows: Data[] }) => {
    const classes = useStyles();
    const [order, setOrder] = useState<Order>('asc');
    const [orderBy, setOrderBy] = useState<keyof Data>('name');
    const [page, setPage] = useState(0);
    const [rowsPerPage, setRowsPerPage] = useState(10);

    const handleRequestSort = (event: React.MouseEvent<unknown>, property: keyof Data) => {
        const isAsc = orderBy === property && order === 'asc';
        batch(() => {
            setOrder(isAsc ? 'desc' : 'asc');
            setOrderBy(property);
        })
    };

    const handleChangePage = (event: unknown, newPage: number) => {
        setPage(newPage);
    };

    const handleChangeRowsPerPage = (event: React.ChangeEvent<HTMLInputElement>) => {
        batch(() => {
            setRowsPerPage(parseInt(event.target.value, 10));
            setPage(0);
        })
    };

    const emptyRows = rowsPerPage - Math.min(rowsPerPage, rows.length - page * rowsPerPage);

    return (
        <div className={classes.root}>
            <TableContainer>
                <Table
                    className={classes.table}
                    aria-labelledby="tableTitle"
                    size={'small'}
                    aria-label="enhanced table"
                >
                    <EnhancedTableHead
                        classes={classes}
                        order={order}
                        orderBy={orderBy}
                        onRequestSort={handleRequestSort}
                    />
                    <TableBody>
                        {stableSort(rows, getComparator(order, orderBy))
                            .slice(page * rowsPerPage, page * rowsPerPage + rowsPerPage)
                            .map((row, index) => {
                                const labelId = `enhanced-table-checkbox-${index}`;

                                return (
                                    <TableRow
                                        tabIndex={-1}
                                        key={row.name}
                                        className={row.implicit ? classes.implicitData : ""}
                                    >
                                        <TableCell component="th" id={labelId} scope="row" padding="none">
                                            {row.name}
                                        </TableCell>
                                        <TableCell align="left">{row.rule}</TableCell>
                                    </TableRow>
                                );
                            })}
                        {emptyRows > 0 && (
                            <TableRow style={{ height: 33 * emptyRows }}>
                                <TableCell colSpan={6} />
                            </TableRow>
                        )}
                    </TableBody>
                </Table>
            </TableContainer>
            <TablePagination
                rowsPerPageOptions={[5, 10, 25]}
                component="div"
                count={rows.length}
                rowsPerPage={rowsPerPage}
                page={page}
                onPageChange={handleChangePage}
                onRowsPerPageChange={handleChangeRowsPerPage}
            />
        </div>
    );
}


export default NamesTable