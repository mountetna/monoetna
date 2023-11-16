import React, { useCallback, useMemo, useEffect, useRef, useState } from 'react';
import { batch } from 'react-redux';
import { ColDef, CellClickedEvent } from 'ag-grid-community';
import { AgGridReact } from 'ag-grid-react';
import { makeStyles } from '@material-ui/core/styles';

import { useWindowDimensions } from '../utils/responsive';

import 'ag-grid-community/styles/ag-grid.css';
import 'ag-grid-community/styles/ag-theme-material.css';



const useStyles = makeStyles((theme) => {

    return {
        container: {
            display: 'flex',
            flexDirection: 'column',
            width: '100%',
            height: '100%',
            border: '1px solid rgb(204, 204, 204)',

            // text input
            '& input[class^=ag-][type=text]:focus': {
                borderBottomColor: 'unset',
            },
            // checked checkbox
            '& .ag-checkbox-input-wrapper.ag-checked::after, & .ag-checkbox-input-wrapper.ag-indeterminate::after, & .ag-radio-button-input-wrapper.ag-checked::after': {
                color: 'rgb(153, 153, 153)',
            },
            // focused cell
            '& .ag-ltr .ag-cell-focus:not(.ag-cell-range-selected):focus-within': {
                borderColor: 'unset',
            },
            // selected row
            '& .ag-row-selected::before, & .ag-row-focus::before, & .ag-row-hover::before': {
                background: 'rgba(153, 153, 153, 0.1)',
            },
        },
        selectedValueContainer: {
            display: 'flex',
            width: '100%',
            alignItems: 'center',
            textAlign: 'left',
            background: 'rgba(153, 153, 153, 0.1)',
        },
        selectedValue: {
            padding: '1em',
        },
        gridContainer: {
            width: '100%',
            height: '100%',
        },
    };
});


const Table = <RowData extends Record<string, any>>({ rows, columns, selectable = false, onCellFocused, onSelectionChanged, className, getRowClass }: {
    rows: RowData[],
    columns: (keyof RowData)[],
    selectable?: boolean,
    onCellFocused?: (data: any) => void,
    onSelectionChanged?: (selection: RowData[]) => void,
    className?: string
    getRowClass?: (data: RowData) => string,
}) => {
    const [focusedCell, setFocusedCell] = useState<string>();

    const classes = useStyles();

    const gridRef = useRef<AgGridReact<RowData>>(null);
    const windowDimensions = useWindowDimensions();

    const commonColumnSettings: ColDef<RowData, any> = {
        filter: true,
        resizable: true,
        headerCheckboxSelection: selectable ? params => {
            const displayedColumns = params.columnApi.getAllDisplayedColumns();
            return displayedColumns[0] === params.column;
        } : false,
        headerCheckboxSelectionFilteredOnly: selectable,
    };
    const [columnDefs, setColumnDefs] = useState<ColDef<RowData, any>[]>(
        // @ts-ignore
        columns.map((col, idx) => ({
            field: col,
            checkboxSelection: selectable && idx == 0,
            ...commonColumnSettings,
        }))
    );

    const defaultColDef = useMemo(() => ({
        sortable: true
    }), []);

    const sizeColumnsToFit = useCallback(() => {
        if (gridRef.current?.api == undefined) { return; }

        gridRef.current.api.sizeColumnsToFit();
    }, []);

    useEffect(sizeColumnsToFit, [windowDimensions]);

    const handleCellClicked = useCallback((event: CellClickedEvent) => {
        const cellValue = event.value;

        batch(() => {
            onCellFocused && onCellFocused(cellValue);
            setFocusedCell(cellValue);
        });
    }, []);

    const handleSelectionChanged = useCallback(() => {
        if (gridRef.current?.api == undefined) { return; }

        onSelectionChanged && onSelectionChanged(
            gridRef.current.api.getSelectedRows()
        );
    }, []);

    return (
        <div className={`${classes.container} ${className ? className : ''}`}>
            <div className={classes.selectedValueContainer}>
                <span className={classes.selectedValue}>
                    {focusedCell ? focusedCell : '‎'}
                </span>
            </div>
            <div className={`ag-theme-material ${classes.gridContainer}`}>
                <AgGridReact
                    ref={gridRef}
                    rowData={rows}
                    columnDefs={columnDefs}
                    defaultColDef={defaultColDef}
                    rowSelection='multiple'
                    suppressRowClickSelection={true}
                    onFirstDataRendered={sizeColumnsToFit}
                    onCellClicked={handleCellClicked}
                    onSelectionChanged={handleSelectionChanged}
                    getRowClass={getRowClass ? params => getRowClass(params.node.data) : undefined}
                />
            </div >
        </div>
    );
};


export default Table;