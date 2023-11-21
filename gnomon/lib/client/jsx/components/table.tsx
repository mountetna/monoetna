import React, { useCallback, useMemo, useEffect, useRef, useState } from 'react';
import { batch } from 'react-redux';
import { ColDef, CellClickedEvent, CellFocusedEvent } from 'ag-grid-community';
import { AgGridReact } from 'ag-grid-react';
import { makeStyles } from '@material-ui/core/styles';
// @ts-ignore
import Color from 'color';

import { useWindowDimensions } from '../utils/responsive';

import 'ag-grid-community/styles/ag-grid.css';
import 'ag-grid-community/styles/ag-theme-material.css';



const useStyles = makeStyles((theme) => {
    const secondaryColor = Color(theme.palette.secondary.main);

    return {
        container: {
            display: 'flex',
            flexDirection: 'column',
            width: '100%',
            height: '100%',
            border: '1px solid #ccc',

            // text input
            '& input[class^=ag-][type=text]:focus': {
                borderBottomColor: 'unset',
            },
            // checked checkbox
            '& .ag-checkbox-input-wrapper.ag-checked::after, & .ag-checkbox-input-wrapper.ag-indeterminate::after, & .ag-radio-button-input-wrapper.ag-checked::after': {
                color: secondaryColor.string(),
            },
            // focused cell
            '& .ag-ltr .ag-cell-focus:not(.ag-cell-range-selected):focus-within': {
                borderColor: 'transparent',
            },
            // selected row
            '& .ag-row-selected::before, & .ag-row-focus::before, & .ag-row-hover::before': {
                background: secondaryColor.alpha(0.1).string(),
            },
            // popups and field pickers
            '& .ag-popup-child:not(.ag-tooltip-custom), & .ag-picker-field-wrapper:focus-within': {
                boxShadow: 'none',
                border: `1px solid ${secondaryColor.string()}`,
            }
        },
        focusedValueContainer: {
            display: 'flex',
            width: '100%',
            alignItems: 'center',
            textAlign: 'left',
            background: secondaryColor.alpha(0.1).string(),
        },
        focusedValue: {
            padding: '1em',
        },
        gridContainer: {
            width: '100%',
            height: '100%',
        },
    };
});


const Table = <RowData extends Record<string, any>>({ rows, columns, selectable = false, showFocusedCell = false, onCellFocused, onSelectionChanged, className, getRowClass }: {
    rows: RowData[],
    columns: (keyof RowData)[],
    selectable?: boolean,
    showFocusedCell?: boolean,
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
        const grid = gridRef.current?.api;
        if (grid == undefined) { return; }

        grid.sizeColumnsToFit();
    }, []);

    useEffect(sizeColumnsToFit, [windowDimensions]);

    const handleCellClicked = useCallback((event: CellClickedEvent | CellFocusedEvent) => {
        const grid = gridRef.current?.api;
        if (grid == undefined) { return; }

        if (!showFocusedCell) {
            grid.clearFocusedCell();
        }

        const cell = grid.getFocusedCell();
        if (cell === null) {
            return;
        }
        const row = grid.getDisplayedRowAtIndex(cell.rowIndex);
        if (row === undefined) {
            return;
        }
        const cellValue = grid.getValue(cell.column, row);

        batch(() => {
            onCellFocused && onCellFocused(cellValue);
            setFocusedCell(cellValue);
        });
    }, []);

    const handleSelectionChanged = useCallback(() => {
        const grid = gridRef.current?.api;
        if (grid == undefined) { return; }

        onSelectionChanged && onSelectionChanged(
            grid.getSelectedRows()
        );
    }, []);

    return (
        <div className={`${classes.container} ${className ? className : ''}`}>
            {showFocusedCell && <div className={classes.focusedValueContainer}>
                <span className={classes.focusedValue}>
                    {focusedCell ? focusedCell : '‎'}
                </span>
            </div>}
            <div className={`ag-theme-material ${classes.gridContainer}`}>
                <AgGridReact
                    ref={gridRef}
                    rowData={rows}
                    columnDefs={columnDefs}
                    defaultColDef={defaultColDef}
                    rowSelection={selectable ? 'multiple' : undefined}
                    suppressRowClickSelection={true}
                    onFirstDataRendered={sizeColumnsToFit}
                    
                    onCellClicked={handleCellClicked}
                    onCellFocused={handleCellClicked}
                    onSelectionChanged={handleSelectionChanged}
                    getRowClass={getRowClass ? params => getRowClass(params.node.data) : undefined}
                />
            </div >
        </div>
    );
};


export default Table;