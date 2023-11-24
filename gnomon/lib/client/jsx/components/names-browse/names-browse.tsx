import React, { useEffect, useState } from 'react';
import { useSelector } from 'react-redux';
import { makeStyles } from '@material-ui/core/styles';

import ProjectHeader from 'etna-js/components/project-header';

import { selectRulesByName } from '../../selectors/rules';
import { fetchNamesWithRuleAndRegexFromMagma } from '../../utils/names';
import { createFnConcurrencyWrapper } from '../../utils/async';
import { useDispatch } from '../../utils/redux';
import { setMagmaNamesListRequest } from '../../actions/names';
import { fetchAndAddRulesFromMagma } from '../../actions/rules';
import { selectMagmaNamesListsByRuleName } from '../../selectors/names';
import NamesBrowseToolbar from './toolbar';
import Counts, { Count } from '../names-create/counts';
import Table from '../table';



const useStyles = makeStyles((theme) => {
    const gridSize = '70';
    const gridSizeSm = '95';

    return {
        projectAndToolbarContainer: {
            display: 'flex',
            borderBottom: '1px solid #eee',
            '& > :first-child, & > :last-child': {
                flex: '1'
            },
            '& > .placeholder': {
                visibility: 'hidden',
            },
        },
        countsList: {
            marginTop: '0.5em',
            textAlign: 'center'
        },
        selectedCounts: {
            display: 'inline-block',
            '& .count, & .separator': {
                display: 'inline-block',
            },
            '& .separator': {
                margin: '0 0.5em'
            },
            '& .count-not-ready': {
                color: 'red',
            },
        },
        tableContainer: {
            display: 'flex',
            justifyContent: 'center',
            margin: '5.5em 0',
        },
        table: {
            width: `${gridSize}vw`,
            height: `${gridSize}vh`,
            [theme.breakpoints.down('sm')]: {
                width: `${gridSizeSm}vw`,
            },
        },
    };
});


export interface Name {
    name: string
    ruleName: string
    author: string
    createdAt: string
}


const NamesBrowse = ({ project_name }: { project_name: string }) => {
    const classes = useStyles();
    const dispatch = useDispatch();

    const [selection, setSelection] = useState<Name[]>([]);
    const [rowData, setRowData] = useState<Name[]>([]);

    const allRules = Object.keys(useSelector(selectRulesByName));
    const magmaNamesListsByRuleName = useSelector(selectMagmaNamesListsByRuleName);

    useEffect(() => {
        async function addRulesFromMagma() {
            dispatch(await fetchAndAddRulesFromMagma(project_name));
        }
        addRulesFromMagma();
    }, []);

    useEffect(() => {
        const fetchNamesForRule = createFnConcurrencyWrapper(fetchNamesWithRuleAndRegexFromMagma, 4);

        async function fetchNamesFromMagma(ruleName: string) {
            try {
                const response = await fetchNamesForRule(project_name, ruleName);

                dispatch(setMagmaNamesListRequest(
                    { status: 'success', response },
                    ruleName,
                ));
            } catch (error) {
                dispatch(setMagmaNamesListRequest(
                    { status: 'error', statusMessage: String(error) },
                    ruleName,
                ));
            }
        }

        allRules.forEach(ruleName => fetchNamesFromMagma(ruleName));
    }, [allRules.length]);

    useEffect(() => {
        const newRowData: Name[] = [];

        for (const [ruleName, namesListRequest] of Object.entries(magmaNamesListsByRuleName)) {
            if (namesListRequest.status != 'success' || namesListRequest.response == undefined) {
                continue;
            }

            for (const magmaName of namesListRequest.response) {

                newRowData.push({
                    name: magmaName.identifier,
                    ruleName,
                    author: magmaName.author,
                    createdAt: magmaName.name_created_at,
                });
            }
        }

        setRowData(newRowData);
    }, [magmaNamesListsByRuleName]);


    const renderCounts = () => {
        const counts: Count[] = [{
            name: 'selected',
            description: 'selected',
            value: selection.length,
            hideAtZero: true,
        }];

        return (
            <div className={classes.countsList}>
                <Counts
                    counts={counts}
                    className={classes.selectedCounts}
                />
            </div>
        );
    };

    return (
        <React.Fragment>
            <div className={classes.projectAndToolbarContainer}>
                <ProjectHeader project_name={project_name} />
                <NamesBrowseToolbar
                    exportData={selection.length ? selection : rowData}
                    exportButtonText={`Export${selection.length ? ' Selection' : ''}`}
                />
                <ProjectHeader project_name={project_name} className="placeholder" />
            </div>
            {renderCounts()}
            <div className={classes.tableContainer}>
                <Table
                    rows={rowData}
                    columns={['name', 'ruleName', 'author', 'createdAt']}
                    selectable={true}
                    showFocusedCell={true}
                    onSelectionChanged={setSelection}
                    className={classes.table}
                    dataTypeLabel='name'
                />
            </div>
        </React.Fragment>
    );
};


export default NamesBrowse;