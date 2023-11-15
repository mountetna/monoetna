import React, { useCallback, useMemo, useEffect, useRef, useState } from 'react'
import { useSelector } from 'react-redux'
import { ColDef, CellClickedEvent, SelectionChangedEvent } from 'ag-grid-community'
import { AgGridReact } from 'ag-grid-react'
import { makeStyles } from '@material-ui/core/styles';

import { selectRulesByName } from '../../selectors/rules'
import { fetchNamesWithRuleAndRegexFromMagma } from '../../utils/names'
import { createFnConcurrencyWrapper } from '../../utils/async'
import { selectPathParts } from '../../selectors/location'
import { useDispatch } from '../../utils/redux'
import { setMagmaNamesListRequest } from '../../actions/names'
import { fetchRulesFromMagma } from '../../actions/rules'
import { selectMagmaNamesListsByRuleName } from '../../selectors/names'

import ProjectHeader from "etna-js/components/project-header";

import NamesBrowseToolbar from './toolbar'

import 'ag-grid-community/styles/ag-grid.css'
import 'ag-grid-community/styles/ag-theme-alpine.css'
import Counts, { Count } from '../names-create/counts';



const useStyles = makeStyles((theme) => {
    const width = "70vw"

    return {
        projectAndToolbarContainer: {
            display: "flex",
            borderBottom: "1px solid #ccc",
            "& > :first-child, & > :last-child": {
                flex: "1"
            },
            "& > .placeholder": {
                visibility: "hidden",
            },
        },
        countsList: {
            marginTop: "0.5em",
            textAlign: "center"
        },
        readyCounts: {
            display: "inline-block",
            '& .count, & .separator': {
                display: "inline-block",
            },
            "& .separator": {
                margin: "0 0.5em"
            },
            "& .count-not-ready": {
                color: "red",
            },
        },
        outerContainer: {
            display: "flex",
            flexDirection: "column",
            alignItems: "center",
            margin: "5.5em 0",
        },
        selectedValue: {
            width,
            textAlign: "left",
        },
        gridContainer: {
            width,
            height: "70vh",
        },
    }
})


interface Name {
    name: string
    ruleName: string
    author: string
    createdAt: string
}


const NamesBrowse = ({ project_name }: { project_name: string }) => {
    const classes = useStyles()

    const [focusedCell, setFocusedCell] = useState<string>()
    const [selection, setSelection] = useState<Name[]>([])

    const commonColumnSettings: ColDef<Name, any> = {
        filter: true,
        resizable: true,
        headerCheckboxSelection: params => {
            const displayedColumns = params.columnApi.getAllDisplayedColumns();
            return displayedColumns[0] === params.column;
        },
        headerCheckboxSelectionFilteredOnly: true,
    }

    const columnDefs: ColDef<Name, any>[] = [
        { field: 'name', checkboxSelection: true, ...commonColumnSettings },
        { field: 'ruleName', ...commonColumnSettings },
        { field: 'author', ...commonColumnSettings },
        { field: 'createdAt', ...commonColumnSettings },
    ]

    const gridRef = useRef<AgGridReact<Name>>(null)
    const dispatch = useDispatch()

    const allRules = Object.keys(useSelector(selectRulesByName))
    const projectName = useSelector(selectPathParts)[0]

    const magmaNamesListsByRuleName = useSelector(selectMagmaNamesListsByRuleName)
    const rowData: Name[] = []

    for (const [ruleName, namesListRequest] of Object.entries(magmaNamesListsByRuleName)) {
        if (namesListRequest.status != "success" || namesListRequest.response == undefined) {
            continue
        }

        for (const magmaName of namesListRequest.response) {
            rowData.push({
                name: magmaName.identifier,
                ruleName,
                author: magmaName.author,
                createdAt: magmaName.name_created_at,
            })
        }
    }


    useEffect(() => {
        async function addRulesFromMagma() {
            dispatch(await fetchRulesFromMagma(projectName))
        }
        addRulesFromMagma()
    }, []);

    useEffect(() => {
        const fetchNamesForRule = createFnConcurrencyWrapper(fetchNamesWithRuleAndRegexFromMagma, 4)

        async function fetchNamesFromMagma(ruleName: string) {
            try {
                const response = await fetchNamesForRule(projectName, ruleName)

                dispatch(setMagmaNamesListRequest(
                    { status: "success", response },
                    ruleName,
                ))
            } catch (error) {
                dispatch(setMagmaNamesListRequest(
                    { status: "error", statusMessage: String(error) },
                    ruleName,
                ))
            }
        }

        allRules.forEach(ruleName => fetchNamesFromMagma(ruleName))
    }, [allRules.length])


    const defaultColDef = useMemo(() => ({
        sortable: true
    }), [])

    const onFirstDataRendered = useCallback((params) => {
        gridRef.current && gridRef.current.api.sizeColumnsToFit();
    }, []);

    const handleCellClicked = useCallback((event: CellClickedEvent) => {
        setFocusedCell(event.value)
    }, [])

    const handleSelectionChanged = useCallback(() => {
        if (gridRef.current == undefined) { return }

        setSelection(
            gridRef.current.api.getSelectedRows()
        )
    }, [])

    const renderCounts = () => {
        const counts: Count[] = [{
            name: "selected",
            description: "selected",
            value: selection.length,
            hideAtZero: true,
        }]

        return (
            <div className={classes.countsList}>
                <Counts
                    counts={counts}
                    className={classes.readyCounts}
                    separator="•"
                />
            </div>
        )
    }

    return (
        <React.Fragment>
            <div className={classes.projectAndToolbarContainer}>
                <ProjectHeader project_name={project_name} />
                <NamesBrowseToolbar
                    exportData={selection.length ? selection : rowData}
                    exportButtonText={`Export${selection.length ? " Selection" : ""}`}
                />
                <ProjectHeader project_name={project_name} className="placeholder" />
            </div>
            {renderCounts()}
            <div className={classes.outerContainer}>
                <div className={classes.selectedValue}>{focusedCell ? focusedCell : " "}</div>
                <div className={`ag-theme-alpine ${classes.gridContainer}`}>
                    <AgGridReact
                        ref={gridRef}
                        rowData={rowData}
                        columnDefs={columnDefs}
                        defaultColDef={defaultColDef}
                        rowSelection='multiple'
                        suppressRowClickSelection={true}
                        onFirstDataRendered={onFirstDataRendered}
                        onCellClicked={handleCellClicked}
                        onSelectionChanged={handleSelectionChanged}
                    />
                </div >
            </div>
        </React.Fragment>
    )
}


export default NamesBrowse


